//! G1/G3/G4: independent P1 induction, source work and material-only changes.
//! These are manufactured coefficients, not a measured metal qualification.
#[path = "winding_ohmic/common.rs"]
mod material;
#[path = "planar_eddy/thermal.rs"]
mod thermal;

use fs_em::planar::{
    PlanarMagnetostatics,
    eddy::{EddyConductor, EddyControls, PlanarEddy},
};
use fs_matdb::MaterialCard;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_rep_mesh::{Metric2, TriComplex2};
use material::{bundle, card, limits, linear, near, point};

fn controls() -> EddyControls {
    EddyControls {
        linear: linear(),
        energy_tolerance_j: 1e-8,
        current_balance_tolerance_a: 1e-9,
    }
}

fn square(
    source: &MaterialCard,
    basis: ElectricalPropertyBasis,
    depth: f64,
) -> PlanarMagnetostatics {
    let mesh = TriComplex2::from_indexed_triangles(
        "eddy-current-four-triangle-mode",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    PlanarMagnetostatics::new(
        mesh,
        vec![bundle(source, 300.0, basis); 4],
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn conductors(
    problem: &PlanarMagnetostatics,
    basis: ElectricalPropertyBasis,
) -> Vec<EddyConductor> {
    problem
        .materials()
        .iter()
        .enumerate()
        .map(|(triangle, material)| EddyConductor {
            triangle,
            region: 7,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                material.clone(),
                basis,
            )
            .unwrap(),
        })
        .collect()
}

fn boundary(offset: f64) -> Vec<(usize, f64)> {
    (0..4).map(|vertex| (vertex, offset)).collect()
}

#[test]
fn isolated_conductor_decay_uses_consistent_mass_and_the_floating_electric_potential() {
    let mut outcomes = Vec::new();
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        for sigma in [18.0, 36.0] {
            let source = card(1.0 / sigma, 0.0, 1.0, basis);
            let depth = 0.3;
            let problem = square(&source, basis, depth);
            // The original magnetostatic solve constructs a nonzero initial field.
            let initial = problem
                .solve(&[12.0; 4], &boundary(0.0), linear(), &mut || false)
                .unwrap();
            near(initial.potential_wb_m()[4], 1.0);
            let state = PlanarEddy::new(
                &initial,
                conductors(&problem, basis),
                &[],
                &boundary(0.0),
                0.0,
                controls(),
                &mut || false,
            )
            .unwrap();
            let dt = 0.1;
            let next = state
                .step_to_candidate(
                    dt,
                    state.conductors().to_vec(),
                    &[0.0; 4],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            // Analytic P1 moments on a unit square: integral(N4)=1/3,
            // integral(N4^2)=1/6. Eliminating the region's floating potential
            // gives C44=sigma*depth*(1/6-1/9), not sigma*depth/6.
            let c = sigma * depth / 18.0;
            let k = 4.0 * depth;
            let amplitude = (c - dt * k / 2.0) / (c + dt * k / 2.0);
            let derivative = (amplitude - 1.0) / dt;
            near(next.state().potential_wb_m()[4], amplitude);
            let step = next.step();
            near(step.stored_before_j, 0.5 * k);
            near(step.stored_after_j, 0.5 * k * amplitude * amplitude);
            near(step.heat_j, dt * c * derivative * derivative);
            near(step.stored_before_j - step.stored_after_j, step.heat_j);
            near(step.supplied_work_j, 0.0);
            near(step.energy_residual_j, 0.0);
            assert_eq!(step.region_currents_a.len(), 1);
            assert_eq!(step.region_currents_a[0].0, 7);
            near(step.region_currents_a[0].1, 0.0);
            for cell in &step.cells {
                let face = problem.mesh().faces()[cell.source.triangle];
                for (local, &vertex) in face.iter().enumerate() {
                    let electric = derivative * if vertex == 4 { -2.0 / 3.0 } else { 1.0 / 3.0 };
                    near(cell.electric_field_v_m[local], electric);
                    near(cell.current_density_a_m2[local], sigma * electric);
                }
                near(cell.volume_m3, depth / 4.0);
                near(cell.heat_j, step.heat_j / 4.0);
                near(cell.power_w, cell.heat_j / dt);
                assert_eq!(cell.source.conductor.source_basis(), basis);
                assert_eq!(
                    cell.source.conductor.material().card_identity(),
                    source.content_hash()
                );
            }
            let field = next.state().field(&mut || false).unwrap();
            near(field.energy_j(), step.stored_after_j);
            for (a, b) in field
                .magnetic_flux_density_t()
                .iter()
                .zip(initial.magnetic_flux_density_t())
            {
                for axis in 0..2 {
                    near(a[axis], amplitude * b[axis]);
                }
            }
            near(state.potential_wb_m()[4], 1.0);
            assert_eq!(state.time_s(), 0.0);
            outcomes.push(amplitude);
        }
    }
    assert!(outcomes[1] > outcomes[0]); // Higher conductivity delays diffusion.
    near(outcomes[0], outcomes[2]);
    near(outcomes[1], outcomes[3]);
}

#[test]
fn retained_induction_has_second_order_midpoint_decay_and_pinned_temperature_refresh() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let source = card(1.0 / 18.0, 0.001, 1.0, basis);
    let problem = square(&source, basis, 0.5);
    let initial = problem
        .solve(&[12.0; 4], &boundary(0.0), linear(), &mut || false)
        .unwrap();
    let mut errors = Vec::new();
    for n in [10, 20, 40] {
        let mut state = PlanarEddy::new(
            &initial,
            conductors(&problem, basis),
            &[],
            &boundary(0.0),
            0.0,
            controls(),
            &mut || false,
        )
        .unwrap();
        let mut heat = 0.0;
        for step in 1..=n {
            let candidate = state
                .step_to_candidate(
                    0.5 * step as f64 / n as f64,
                    state.conductors().to_vec(),
                    &[0.0; 4],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            heat += candidate.step().heat_j;
            state = candidate.into_parts().0;
        }
        near(state.energy_j() + heat, initial.energy_j());
        errors.push((state.potential_wb_m()[4] - (-2.0_f64).exp()).abs());
        assert_eq!(state.time_s(), 0.5);
    }
    assert!(errors[0] > 3.8 * errors[1]);
    assert!(errors[1] > 3.8 * errors[2]);

    let state = PlanarEddy::new(
        &initial,
        conductors(&problem, basis),
        &[],
        &boundary(0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let refreshed: Vec<_> = state
        .conductors()
        .iter()
        .map(|old| EddyConductor {
            triangle: old.triangle,
            region: old.region,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                old.conductor
                    .material()
                    .at_query_point(&source, &point(400.0))
                    .unwrap(),
                basis,
            )
            .unwrap(),
        })
        .collect();
    let hot = state
        .step_to_candidate(0.1, refreshed, &[0.0; 4], controls(), &mut || false)
        .unwrap();
    let c = (18.0 / 1.1) * 0.5 / 18.0;
    let k = 4.0 * 0.5;
    near(
        hot.state().potential_wb_m()[4],
        (c - 0.1 * k / 2.0) / (c + 0.1 * k / 2.0),
    );
    near(hot.step().stored_before_j, initial.energy_j());
    for cell in hot.state().conductors() {
        near(cell.conductor.resistivity_ohm_m(), 1.1 / 18.0);
        assert_eq!(
            cell.conductor.material().card_identity(),
            source.content_hash()
        );
    }
}

#[test]
fn balanced_impressed_current_builds_the_same_retained_field_and_then_decays_as_heat() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [
        card(1.0 / 18.0, 0.0, 1.0, basis),
        card(1.0 / 18.0, 0.0, 1.0, basis),
    ];
    let (problem, profiles, boundary_nodes) = material::problem(&cards, basis, 0.3);
    let bc: Vec<_> = boundary_nodes.iter().map(|&node| (node, 0.0)).collect();
    let zero = vec![0.0; problem.mesh().faces().len()];
    let initial = problem.solve(&zero, &bc, linear(), &mut || false).unwrap();
    let state = PlanarEddy::new(
        &initial,
        conductors(&problem, basis),
        &[],
        &bc,
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let candidate = state
        .step_to_candidate(
            0.001,
            state.conductors().to_vec(),
            &profiles[0],
            controls(),
            &mut || false,
        )
        .unwrap();
    let driven = candidate.step();
    assert!(driven.supplied_work_j > 0.0);
    assert!(driven.stored_after_j > 0.0);
    assert!(driven.heat_j > 0.0);
    near(
        driven.supplied_work_j,
        driven.stored_after_j + driven.heat_j,
    );
    let independent_work: f64 = problem
        .mesh()
        .faces()
        .iter()
        .enumerate()
        .map(|(i, face)| {
            // This fixture's uniform 8x8 mesh has triangle volume depth/128.
            // The accepted source starts at zero and ramps to this endpoint.
            0.5 * profiles[0][i] * 0.3 / 128.0
                * face
                    .iter()
                    .map(|&v| candidate.state().potential_wb_m()[v as usize] / 3.0)
                    .sum::<f64>()
        })
        .sum();
    near(driven.supplied_work_j, independent_work);
    let (running, _) = candidate.into_parts();
    let ramp_down = running
        .step_to_candidate(
            0.002,
            running.conductors().to_vec(),
            &zero,
            controls(),
            &mut || false,
        )
        .unwrap();
    near(
        ramp_down.step().stored_after_j - ramp_down.step().stored_before_j
            + ramp_down.step().heat_j,
        ramp_down.step().supplied_work_j,
    );
    assert_eq!(ramp_down.state().impressed_current_density_a_m2(), zero);
    let (unforced, _) = ramp_down.into_parts();
    let coast = unforced
        .step_to_candidate(
            0.003,
            unforced.conductors().to_vec(),
            &zero,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(coast.step().stored_after_j < coast.step().stored_before_j);
    near(
        coast.step().stored_before_j - coast.step().stored_after_j,
        coast.step().heat_j,
    );
    near(coast.step().supplied_work_j, 0.0);
    near(state.energy_j(), 0.0);
}

#[test]
fn insulating_algebraic_field_follows_both_source_endpoints_without_alternating_history() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [
        card(1.0 / 18.0, 0.0, 1.0, basis),
        card(1.0 / 18.0, 0.0, 1.0, basis),
    ];
    let (problem, profiles, nodes) = material::problem(&cards, basis, 0.3);
    let bc: Vec<_> = nodes.into_iter().map(|v| (v, 0.0)).collect();
    let insulators: Vec<_> = (0..problem.mesh().faces().len()).collect();
    let zero = vec![0.0; insulators.len()];
    let initial = problem.solve(&zero, &bc, linear(), &mut || false).unwrap();
    let excited = problem
        .solve(&profiles[0], &bc, linear(), &mut || false)
        .unwrap();
    let mut limits = controls();
    limits.current_balance_tolerance_a = 1e-8;
    // Electrical roles are explicit: this section has no modeled conductor.
    // Every free magnetic coordinate is therefore an algebraic air coordinate.
    assert!(
        PlanarEddy::new(&excited, vec![], &insulators, &bc, 0.0, limits, &mut || {
            false
        })
        .is_err()
    );
    let consistent = PlanarEddy::new_with_impressed_current(
        &excited,
        vec![],
        &insulators,
        &bc,
        &profiles[0],
        0.0,
        limits,
        &mut || false,
    )
    .unwrap();
    assert_eq!(consistent.potential_wb_m(), excited.potential_wb_m());
    let mut state = PlanarEddy::new(&initial, vec![], &insulators, &bc, 0.0, limits, &mut || {
        false
    })
    .unwrap();
    for (time, input, expected) in [
        (0.001, profiles[0].as_slice(), excited.potential_wb_m()),
        (0.002, zero.as_slice(), initial.potential_wb_m()),
    ] {
        let candidate = state
            .step_to_candidate(time, vec![], input, limits, &mut || false)
            .unwrap();
        for (&actual, &reference) in candidate.state().potential_wb_m().iter().zip(expected) {
            near(actual, reference);
        }
        assert_eq!(candidate.state().impressed_current_density_a_m2(), input);
        assert_eq!(candidate.step().heat_j, 0.0);
        near(
            candidate.step().stored_after_j - candidate.step().stored_before_j,
            candidate.step().supplied_work_j,
        );
        state = candidate.into_parts().0;
    }
    near(state.energy_j(), 0.0);
}

#[test]
fn field_and_material_history_survive_source_budget_cancel_and_topology_refusals() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let source = card(1.0 / 18.0, 0.001, 1.0, basis);
    let problem = square(&source, basis, 0.5);
    let initial = problem
        .solve(&[12.0; 4], &boundary(0.0), linear(), &mut || false)
        .unwrap();
    let state = PlanarEddy::new(
        &initial,
        conductors(&problem, basis),
        &[],
        &boundary(0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let old = state.potential_wb_m().to_vec();
    let energy = state.energy_j();
    assert!(
        state
            .step_to_candidate(
                0.1,
                state.conductors().to_vec(),
                &[1.0; 4],
                controls(),
                &mut || false
            )
            .is_err()
    ); // No hidden impressed-current return.
    let mut starved = controls();
    starved.linear.max_iterations = 0;
    assert!(
        state
            .step_to_candidate(
                0.1,
                state.conductors().to_vec(),
                &[0.0; 4],
                starved,
                &mut || false
            )
            .is_err()
    );
    let mut checks = 0;
    assert!(
        state
            .step_to_candidate(
                0.1,
                state.conductors().to_vec(),
                &[0.0; 4],
                controls(),
                &mut || {
                    checks += 1;
                    checks > 8
                }
            )
            .is_err()
    );
    assert!(checks > 8);
    let other = card(1.0 / 36.0, 0.001, 1.0, basis);
    let mut replaced = state.conductors().to_vec();
    replaced[0].conductor =
        ResolvedConductorState::try_from_material_with_basis(bundle(&other, 300.0, basis), basis)
            .unwrap();
    assert!(
        state
            .step_to_candidate(0.1, replaced, &[0.0; 4], controls(), &mut || false)
            .is_err()
    );
    let mut disconnected = state.conductors().to_vec();
    disconnected[1].region = 8;
    disconnected[3].region = 8;
    assert!(
        PlanarEddy::new(
            &initial,
            disconnected,
            &[],
            &boundary(0.0),
            0.0,
            controls(),
            &mut || false
        )
        .is_err()
    ); // Opposite triangles share only a point.
    assert_eq!(state.potential_wb_m(), old);
    near(state.energy_j(), energy);
    assert_eq!(state.time_s(), 0.0);
    let a = state
        .step_to_candidate(
            0.1,
            state.conductors().to_vec(),
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    let b = state
        .step_to_candidate(
            0.1,
            state.conductors().to_vec(),
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(a.state().potential_wb_m(), b.state().potential_wb_m());
    near(a.step().heat_j, b.step().heat_j);
}

