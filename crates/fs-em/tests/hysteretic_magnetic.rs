//! G1/G3/G4: actual planar field histories, return-point memory and work balance.
//! Synthetic play coefficients verify the implementation, not catalog materials.

#![cfg(feature = "magnetic-hysteresis")]

#[path = "hysteretic_magnetic/path.rs"]
mod path;

use fs_em::planar::{
    PlanarMagnetostatics,
    nonlinear::{
        NonlinearMagneticControls,
        hysteretic::{PlanarHystereticControls, PlanarHystereticMagnetostatics},
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        hysteresis::{
            source::{
                HysteresisBranchSources, HysteresisSources, ResolvedScalarHysteresis,
                ResolvedVectorHysteresis,
            },
            vector::VectorHysteresisControls,
        },
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;
use std::sync::Arc;

fn card(values: Vec<(PropertyKey, f64)>) -> Arc<MaterialCard> {
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
                    source: "manufactured isotropic magnetic play fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic magnetic law".into(),
                phase: "solid".into(),
                process: "manufactured field history".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}

fn source(mu: f64, pinning: bool) -> ResolvedVectorHysteresis {
    let mut roles = HysteresisSources {
        background_permeability: PropertyKey::new("background_mu", MAGNETIC_PERMEABILITY_DIMS),
        field_limit: PropertyKey::new("field_limit", FIELD_STRENGTH_DIMS),
        branches: vec![],
    };
    let mut values = vec![
        (roles.background_permeability.clone(), mu),
        (roles.field_limit.clone(), 30.0),
    ];
    if pinning {
        let branch = HysteresisBranchSources {
            susceptibility: PropertyKey::new("branch_slope", MAGNETIC_PERMEABILITY_DIMS),
            pinning_field: PropertyKey::new("branch_pinning", FIELD_STRENGTH_DIMS),
            saturation_polarization: PropertyKey::new("branch_saturation", FLUX_DENSITY_DIMS),
        };
        values.extend([
            (branch.susceptibility.clone(), 0.2),
            (branch.pinning_field.clone(), 2.0),
            (branch.saturation_polarization.clone(), 1.0),
        ]);
        roles.branches.push(branch);
    }
    let scalar = ResolvedScalarHysteresis::resolve(
        card(values),
        &QueryPoint::new(),
        roles,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedVectorHysteresis::from_isotropic_source(scalar).unwrap()
}

fn linear_material(mu: f64) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        &card(vec![(
            PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
            mu,
        )]),
        &QueryPoint::new(),
        &[ScalarPropertyRequirement::try_new(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn vertices() -> Vec<[f64; 2]> {
    vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]]
}

fn square(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "magnetic-history-square",
        vertices(),
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}

fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 100,
        max_dofs: 100,
        max_triplets: 10000,
    }
}

fn controls() -> PlanarHystereticControls {
    PlanarHystereticControls {
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
    }
}

fn boundary(gauge: f64) -> Vec<(usize, f64)> {
    (0..4).map(|i| (i, gauge)).collect()
}

fn problem(depth: f64) -> PlanarHystereticMagnetostatics {
    PlanarHystereticMagnetostatics::new(
        square(depth),
        vec![source(0.01, true); 4],
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-8 * (1.0 + b.abs()), "{a:.16e} != {b:.16e}");
}

#[test]
fn field_cycles_follow_analytic_play_and_retain_nontrivial_minor_loops() {
    let depth = 0.7;
    let gauge = 1.25;
    let owner = problem(depth);
    let mut state = owner.demagnetized(gauge, &mut || false).unwrap();
    let directions = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]];
    let mut old_z: f64 = 0.0;
    let mut old_b = 0.0;
    let mut stored = 0.0;
    let mut minor_return = None;
    // After saturation, 0 -> 6 -> 1 -> 5.5 -> 1 -> 6 -> 0 contains two
    // nontrivial nested play loops and returns to the same remanent state.
    for (index, h) in [12.0, 0.0, 6.0, 1.0, 5.5, 1.0, 6.0, 0.0, -12.0, 0.0]
        .into_iter()
        .enumerate()
    {
        // Four equal triangles: integral H.curl(N_center)=2*depth*H,
        // integral J*N_center=depth*J/3. Hence H=J/6 independently of mu.
        let z = old_z.clamp(0.2 * (h - 2.0), 0.2 * (h + 2.0)).clamp(-1.0, 1.0);
        let b = 0.01 * h + z;
        let energy = depth * (0.005 * h * h + z * z / 0.4);
        let candidate = owner
            .candidate(&state, &[6.0 * h; 4], &boundary(gauge), controls(), &mut || false)
            .unwrap();
        let next = candidate.state();
        let receipt = candidate.step();
        near(next.potential_wb_m()[4] - gauge, 0.5 * b);
        near(next.stored_energy_j(), energy);
        near(receipt.stored_energy_change_j, energy - stored);
        near(receipt.irreversible_heat_j, depth * 2.0 * (z - old_z).abs());
        near(receipt.supplied_work_j, depth * h * (b - old_b));
        near(receipt.impressed_current_work_j, receipt.supplied_work_j);
        near(receipt.boundary_work_j, 0.0);
        near(receipt.field_work_j, receipt.supplied_work_j);
        near(
            receipt.supplied_work_j,
            receipt.stored_energy_change_j
                + receipt.irreversible_heat_j
                + receipt.numerical_dissipation_j,
        );
        near(receipt.energy_residual_j, 0.0);
        assert!(receipt.numerical_dissipation_j >= 0.0);
        assert_eq!(receipt.cells.len(), 4);
        for (cell, direction) in directions.iter().enumerate() {
            let physical = &next.cell_states()[cell];
            for axis in 0..2 {
                near(next.magnetic_flux_density_t()[cell][axis], b * direction[axis]);
                near(physical.field_a_m()[axis], h * direction[axis]);
                near(physical.polarizations_t()[0][axis], z * direction[axis]);
            }
            near(physical.field_a_m()[2], 0.0);
            near(receipt.cells[cell].volume_m3, depth / 4.0);
            near(receipt.cells[cell].irreversible_heat_j, depth * 0.5 * (z - old_z).abs());
            assert_eq!(receipt.cells[cell].triangle, cell);
            assert!(receipt.cells[cell].flux_residual_t <= 1e-10);
        }
        if index == 1 {
            assert!(next.magnetic_flux_density_t()[0][0] > 0.39); // Remanence at J=0.
            minor_return = Some(next.magnetic_flux_density_t().to_vec());
        }
        if index == 7 {
            for (actual, returned) in next
                .magnetic_flux_density_t()
                .iter()
                .zip(minor_return.as_ref().unwrap())
            {
                near(actual[0], returned[0]);
                near(actual[1], returned[1]);
            }
        }
        (state, _) = candidate.into_parts();
        old_z = z;
        old_b = b;
        stored = energy;
    }
}

fn full_loop(segments: usize) -> (f64, f64, f64) {
    let owner = problem(0.7);
    let virgin = owner.demagnetized(0.0, &mut || false).unwrap();
    let (mut state, _) = owner
        .candidate(&virgin, &[72.0; 4], &boundary(0.0), controls(), &mut || false)
        .unwrap()
        .into_parts();
    let start_energy = state.stored_energy_j();
    let (mut heat, mut numerical, mut supplied) = (0.0, 0.0, 0.0);
    for (from, to) in [(12.0, -12.0), (-12.0, 12.0)] {
        for i in 1..=segments {
            let h = from + (to - from) * i as f64 / segments as f64;
            let (next, receipt) = owner
                .candidate(&state, &[6.0 * h; 4], &boundary(0.0), controls(), &mut || false)
                .unwrap()
                .into_parts();
            heat += receipt.irreversible_heat_j;
            numerical += receipt.numerical_dissipation_j;
            supplied += receipt.supplied_work_j;
            near(receipt.energy_residual_j, 0.0);
            state = next;
        }
    }
    near(state.stored_energy_j(), start_energy);
    near(supplied, heat + numerical);
    (heat, numerical, supplied)
}

#[test]
fn loop_refinement_preserves_physical_heat_and_removes_endpoint_numerical_loss() {
    let coarse = full_loop(8);
    let fine = full_loop(64);
    // Once saturated, total polarization variation is 4*zs per complete loop.
    let exact_heat = 0.7 * 4.0 * 2.0 * 1.0;
    near(coarse.0, exact_heat);
    near(fine.0, exact_heat);
    assert!(fine.1 > 0.0 && fine.1 < 0.2 * coarse.1);
    assert!(fine.2 < coarse.2 && fine.2 > exact_heat);
}

#[test]
fn branchless_heterogeneous_materials_match_existing_linear_field_solver() {
    let permeability = [0.4, 0.7, 1.2, 2.0];
    let reference = PlanarMagnetostatics::new(
        square(0.7),
        permeability.into_iter().map(linear_material).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let owner = PlanarHystereticMagnetostatics::new(
        square(0.7),
        permeability.into_iter().map(|mu| source(mu, false)).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let gauge = 2.0;
    let mut state = owner.demagnetized(gauge, &mut || false).unwrap();
    for current in [12.0, -6.0, 0.0, 12.0] {
        let exact = reference
            .solve(&[current; 4], &boundary(gauge), controls().nonlinear.linear, &mut || false)
            .unwrap();
        let (next, receipt) = owner
            .candidate(&state, &[current; 4], &boundary(gauge), controls(), &mut || false)
            .unwrap()
            .into_parts();
        for (&a, &b) in next.potential_wb_m().iter().zip(exact.potential_wb_m()) {
            near(a, b);
        }
        for (a, b) in next.magnetic_flux_density_t().iter().zip(exact.magnetic_flux_density_t()) {
            near(a[0], b[0]);
            near(a[1], b[1]);
        }
        near(next.stored_energy_j(), exact.energy_j());
        assert_eq!(receipt.irreversible_heat_j, 0.0);
        near(receipt.energy_residual_j, 0.0);
        state = next;
    }
}

#[test]
fn rotating_boundary_field_uses_native_vector_memory_and_actual_reaction_work() {
    let depth = 0.4;
    let material = source(0.01, true);
    // Every vertex is on the exterior trace; an interior degree of freedom
    // cannot be prescribed through the original physical boundary API.
    let mesh = TriComplex2::from_indexed_triangles(
        "rotating-boundary-hysteresis",
        vertices()[..4].to_vec(),
        vec![[0, 1, 2], [0, 2, 3]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    let owner = PlanarHystereticMagnetostatics::new(
        mesh,
        vec![material.clone(); 2],
        limits(),
        &mut || false,
    )
    .unwrap();
    let gauge = 0.75;
    let mut state = owner.demagnetized(gauge, &mut || false).unwrap();
    let mut memory = material.law().initial_state();
    for b in [
        [0.8, 0.0], [0.6, 0.8], [-0.6, 0.8], [-0.8, 0.0], [0.0, -0.8], [0.8, 0.0],
    ] {
        let prescribed: Vec<_> = vertices()
            .iter()
            .take(4)
            .enumerate()
            .map(|(i, &[x, y])| (i, gauge + b[0] * y - b[1] * x))
            .collect();
        let native = material
            .law()
            .candidate_flux_density(
                &memory, [b[0], b[1], 0.0], controls().material, &mut || false,
            )
            .unwrap();
        let (next, receipt) = owner
            .candidate(&state, &[0.0; 2], &prescribed, controls(), &mut || false)
            .unwrap()
            .into_parts();
        for cell in next.cell_states() {
            for axis in 0..3 {
                near(cell.field_a_m()[axis], native.state.field_a_m()[axis]);
                near(cell.polarizations_t()[0][axis], native.state.polarizations_t()[0][axis]);
            }
        }
        near(next.stored_energy_j(), depth * native.state.stored_energy_j_m3());
        near(receipt.irreversible_heat_j, depth * native.irreversible_heat_j_m3);
        near(receipt.numerical_dissipation_j, depth * native.numerical_dissipation_j_m3);
        // With J=0 and every A prescribed, all supply comes from boundary reactions.
        near(receipt.supplied_work_j, depth * native.endpoint_work_j_m3);
        near(receipt.impressed_current_work_j, 0.0);
        near(receipt.boundary_work_j, receipt.supplied_work_j);
        near(receipt.energy_residual_j, 0.0);
        assert_eq!(receipt.nonlinear.newton_iterations, 0);
        memory = native.state;
        state = next;
    }
}

#[test]
fn failed_and_cancelled_field_trials_never_commit_material_history() {
    let owner = problem(1.0);
    let virgin = owner.demagnetized(0.0, &mut || false).unwrap();
    let (state, _) = owner
        .candidate(&virgin, &[24.0; 4], &boundary(0.0), controls(), &mut || false)
        .unwrap()
        .into_parts();
    let memory = state.cell_states().to_vec();
    let potential = state.potential_wb_m().to_vec();
    let stored = state.stored_energy_j().to_bits();
    let mut exhausted = controls();
    exhausted.nonlinear.max_newton_iterations = 0;
    assert!(
        owner.candidate(&state, &[-72.0; 4], &boundary(0.0), exhausted, &mut || false)
            .is_err()
    );
    assert!(
        owner.candidate(&state, &[24.0; 4], &boundary(0.0), controls(), &mut || true)
            .is_err()
    );
    assert!(
        owner.candidate(&state, &[f64::NAN; 4], &boundary(0.0), controls(), &mut || false)
            .is_err()
    );
    assert!(
        owner.candidate(&state, &[240.0; 4], &boundary(0.0), controls(), &mut || false)
            .is_err()
    );
    let mut calls = 0;
    let accepted = owner
        .candidate(&state, &[-30.0; 4], &boundary(0.0), controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 10);
    let mut counter = 0;
    assert!(
        owner.candidate(&state, &[-30.0; 4], &boundary(0.0), controls(), &mut || {
            counter += 1;
            counter == calls / 2
        }).is_err()
    );
    assert_eq!(state.cell_states(), memory);
    assert_eq!(state.potential_wb_m(), potential);
    assert_eq!(state.stored_energy_j().to_bits(), stored);
    let repeated = owner
        .candidate(&state, &[-30.0; 4], &boundary(0.0), controls(), &mut || false)
        .unwrap();
    assert_eq!(accepted.state().cell_states(), repeated.state().cell_states());
    assert_eq!(accepted.state().potential_wb_m(), repeated.state().potential_wb_m());
    assert!(
        PlanarHystereticMagnetostatics::new(
            square(1.0), vec![source(0.01, true); 3], limits(), &mut || false,
        ).is_err()
    );
}

#[test]
fn explicit_initial_guess_admits_valid_boundary_loading_near_source_limit() {
    let owner = problem(0.6);
    let previous = owner.demagnetized(0.0, &mut || false).unwrap();
    let initial: Vec<_> = vertices().iter().map(|p| 1.25 * p[1]).collect();
    let prescribed: Vec<_> = (0..4).map(|i| (i, initial[i])).collect();
    // The unchanged interior A=0 with new top-boundary A=1.25 would give B=2.5
    // on one cell, outside the fitted maximum 0.01*30+1=1.3 T. The affine
    // solution has B=1.25 and H=25 A/m throughout, within the original domain.
    assert!(
        owner.candidate(&previous, &[0.0; 4], &prescribed, controls(), &mut || false)
            .is_err()
    );
    let candidate = owner
        .candidate_from_initial(
            &previous,
            &[0.0; 4],
            &prescribed,
            &initial,
            controls(),
            &mut || false,
        )
        .unwrap();
    for cell in candidate.state().cell_states() {
        near(cell.flux_density_t()[0], 1.25);
        near(cell.flux_density_t()[1], 0.0);
        near(cell.field_a_m()[0], 25.0);
        near(cell.polarizations_t()[0][0], 1.0);
    }
    near(candidate.step().irreversible_heat_j, 0.6 * 2.0);
    near(candidate.step().boundary_work_j, 0.6 * 25.0 * 1.25);
    near(candidate.step().energy_residual_j, 0.0);
    assert_eq!(candidate.step().nonlinear.newton_iterations, 0);
    assert_eq!(previous.stored_energy_j(), 0.0);
}

#[test]
fn disconnected_gauge_changes_preserve_each_history_and_do_no_physical_work() {
    let mut coordinates = vertices();
    coordinates.extend(vertices().iter().map(|&[x, y]| [x + 2.0, y]));
    let mut triangles = vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
    triangles.extend([[5, 6, 9], [6, 7, 9], [7, 8, 9], [8, 5, 9]]);
    let mesh = TriComplex2::from_indexed_triangles(
        "disconnected-hysteretic-fields",
        coordinates,
        triangles,
        Metric2::planar(0.7).unwrap(),
    )
    .unwrap();
    let owner = PlanarHystereticMagnetostatics::new(
        mesh, vec![source(0.01, true); 8], limits(), &mut || false,
    )
    .unwrap();
    let currents = [18.0, 18.0, 18.0, 18.0, -24.0, -24.0, -24.0, -24.0];
    let boundary_vertices = [0, 1, 2, 3, 5, 6, 7, 8];
    let zero_boundary: Vec<_> = boundary_vertices.map(|i| (i, 0.0)).into();
    let virgin = owner.demagnetized(0.0, &mut || false).unwrap();
    let (old, _) = owner
        .candidate(&virgin, &currents, &zero_boundary, controls(), &mut || false)
        .unwrap()
        .into_parts();
    let shifted: Vec<_> = boundary_vertices
        .map(|i| (i, if i < 5 { 7.0 } else { -11.0 }))
        .into();
    let candidate = owner
        .candidate(&old, &currents, &shifted, controls(), &mut || false)
        .unwrap();
    for (i, (&a, &b)) in candidate
        .state()
        .potential_wb_m()
        .iter()
        .zip(old.potential_wb_m())
        .enumerate()
    {
        near(a - b, if i < 5 { 7.0 } else { -11.0 });
    }
    for (after, before) in candidate
        .state()
        .cell_states()
        .iter()
        .zip(old.cell_states())
    {
        for axis in 0..3 {
            near(after.flux_density_t()[axis], before.flux_density_t()[axis]);
            near(after.polarizations_t()[0][axis], before.polarizations_t()[0][axis]);
        }
    }
    near(candidate.state().stored_energy_j(), old.stored_energy_j());
    // Both reported pairings use the component's own gauge reference.
    near(candidate.step().impressed_current_work_j, 0.0);
    near(candidate.step().boundary_work_j, 0.0);
    near(candidate.step().supplied_work_j, 0.0);
    near(candidate.step().irreversible_heat_j, 0.0);
    near(candidate.step().numerical_dissipation_j, 0.0);
    near(candidate.step().energy_residual_j, 0.0);
}
