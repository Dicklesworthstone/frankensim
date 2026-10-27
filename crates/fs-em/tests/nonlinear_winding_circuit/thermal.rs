//! G1/G3/G4: source-pinned mu(T), native voltage drive and signed magnetic work.
//! These manufactured laws verify the discrete equations, not measured metals.
use super::common;
#[path = "../winding_ohmic/common.rs"]
mod source;

use common::{controls, near};
use fs_em::{
    planar::{
        eddy::{EddyConductor, nonlinear::PlanarNonlinearEddy},
        nonlinear::{IsotropicMagneticMaterial, PlanarNonlinearMagnetostatics},
    },
    winding::{
        nonlinear::{
            NonlinearWindingCircuit, NonlinearWindingError, NonlinearWindingInitial,
            NonlinearWindings,
        },
        ohmic::WindingConductor,
    },
};
use fs_matdb::MaterialCard;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_phs::circuit::{Branch, CircuitGraph};

fn materials(card: &MaterialCard, t: f64, count: usize) -> Vec<IsotropicMagneticMaterial> {
    vec![
        IsotropicMagneticMaterial::Linear(source::bundle(
            card,
            t,
            ElectricalPropertyBasis::Resistivity,
        ));
        count
    ]
}
fn conductors(materials: &[IsotropicMagneticMaterial]) -> Vec<WindingConductor> {
    materials
        .iter()
        .enumerate()
        .map(|(triangle, material)| WindingConductor {
            triangle,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                material.material().clone(),
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap(),
        })
        .collect()
}
fn initial(problem: &PlanarNonlinearMagnetostatics) -> NonlinearWindingCircuit<'_> {
    let profiles = common::profiles(1);
    NonlinearWindingCircuit::new(
        common::initial_field(problem, &profiles, &[1.0]),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 999.0 }),
                (2, 0, Branch::Resistor { ohms: 0.3 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.0],
        },
        conductors(problem.materials()),
        NonlinearWindingInitial {
            currents_a: vec![1.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![1.5],
        },
        common::config(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn permeability_temperature_changes_native_current_field_and_signed_work_without_reset() {
    let depth = 0.3;
    let dt = 0.2;
    let voltage = 9.0;
    let resistance = 72.0 * depth / 18.0 + 0.3;
    for slope in [-0.0002, 0.0, 0.001] {
        let card = source::affine_magnetic_card(
            1.0 / 18.0,
            0.0,
            1.0,
            slope,
            ElectricalPropertyBasis::Resistivity,
        );
        let problem = PlanarNonlinearMagnetostatics::new(
            common::mesh(depth),
            materials(&card, 300.0, 8),
            source::limits(),
            &mut || false,
        )
        .unwrap();
        let original = initial(&problem);
        let at = materials(&card, 500.0, 8);
        let next = original
            .step_with_materials_candidate(
                dt,
                at.clone(),
                vec![],
                conductors(&at),
                &[voltage],
                source::limits(),
                controls(),
                &mut || false,
            )
            .unwrap();
        let mu = 1.0 + 200.0 * slope;
        // Independent RL reduction: lambda0=2*depth and lambda1=2*depth*mu*I1.
        let current = (voltage + 2.0 * depth / dt) / (resistance + 2.0 * depth * mu / dt);
        let b = mu * current;
        let parameter = depth * (1.0 / mu - 1.0);
        let numerical = depth * (b - 1.0).powi(2) / mu;
        near(next.circuit().currents_a()[0], current);
        near(next.circuit().flux_linkages_wb()[0], 2.0 * depth * b);
        near(next.circuit().field().potential_wb_m()[6], b / 2.0);
        near(next.circuit().field().potential_wb_m()[7], -b / 2.0);
        let step = next.step();
        near(step.stored_before_j, depth);
        near(step.stored_after_j, depth * b * b / mu);
        near(step.parameter_work_j, parameter);
        near(step.field.parameter_work_j, parameter);
        near(step.numerical_dissipation_j, numerical);
        near(step.heat_j, dt * resistance * current * current);
        near(step.supplied_work_j, dt * voltage * current);
        near(step.energy_residual_j, 0.0);
        near(
            step.field
                .magnetic_cells
                .iter()
                .map(|c| c.parameter_work_j)
                .sum(),
            parameter,
        );
        for cell in &step.field.magnetic_cells {
            near(cell.parameter_work_j, parameter / 8.0);
        }
        near(
            step.stored_after_j - step.stored_before_j + step.heat_j + numerical,
            step.supplied_work_j + parameter,
        );

        // Owned refreshed geometry survives clone and decomposition. The next
        // interval uses its accepted lambda, not a newly initialized current.
        let (retained, _) = next.into_parts();
        let restarted = retained.clone();
        let at = materials(&card, 600.0, 8);
        let second = restarted
            .step_with_materials_candidate(
                2.0 * dt,
                at.clone(),
                vec![],
                conductors(&at),
                &[voltage],
                source::limits(),
                controls(),
                &mut || false,
            )
            .unwrap();
        let mu2 = 1.0 + 300.0 * slope;
        let current2 = (voltage + 2.0 * depth * b / dt) / (resistance + 2.0 * depth * mu2 / dt);
        near(second.circuit().currents_a()[0], current2);
        near(
            second.step().parameter_work_j,
            depth * b * b * (1.0 / mu2 - 1.0 / mu),
        );
        near(second.step().stored_before_j, retained.field().energy_j());
        near(second.step().energy_residual_j, 0.0);
        near(original.time_s(), 0.0);
        near(original.currents_a()[0], 1.0);
        near(original.field().potential_wb_m()[6], 0.5);
    }
}

#[test]
fn conducting_core_refresh_closes_local_joule_numerical_and_parameter_accounts() {
    let depth = 0.3;
    let dt = 0.1;
    let card = source::affine_magnetic_card(
        1.0 / 18.0,
        0.0,
        1.0,
        0.001,
        ElectricalPropertyBasis::Resistivity,
    );
    let problem = PlanarNonlinearMagnetostatics::new(
        common::material::mesh(depth),
        materials(&card, 300.0, 4),
        source::limits(),
        &mut || false,
    )
    .unwrap();
    let boundary = common::material::boundary(0.0);
    let start = problem
        .solve(
            &[12.0; 4],
            &boundary,
            &[0.0; 5],
            controls().field.nonlinear,
            &mut || false,
        )
        .unwrap();
    let eddy_sources = |at: &[IsotropicMagneticMaterial]| {
        conductors(at)
            .into_iter()
            .map(|c| EddyConductor {
                triangle: c.triangle,
                region: 7,
                conductor: c.conductor,
            })
            .collect()
    };
    let original = PlanarNonlinearEddy::new(
        &start,
        eddy_sources(problem.materials()),
        &[],
        &boundary,
        0.0,
        controls().field,
        &mut || false,
    )
    .unwrap();
    let at = materials(&card, 500.0, 4);
    let next = original
        .step_with_materials_candidate(
            dt,
            at.clone(),
            eddy_sources(&at),
            &[0.0; 4],
            source::limits(),
            controls().field,
            &mut || false,
        )
        .unwrap();
    let k0 = 4.0 * depth;
    let k1 = k0 / 1.2;
    let c = depth;
    let a = c / (c + dt * k1);
    near(next.state().potential_wb_m()[4], a);
    near(next.step().stored_before_j, 0.5 * k0);
    near(next.step().stored_after_j, 0.5 * k1 * a * a);
    near(next.step().heat_j, c * (a - 1.0).powi(2) / dt);
    near(
        next.step().numerical_dissipation_j,
        0.5 * k1 * (a - 1.0).powi(2),
    );
    near(next.step().parameter_work_j, 0.5 * (k1 - k0));
    near(next.step().supplied_work_j, 0.0);
    near(next.step().energy_residual_j, 0.0);
    for cell in &next.step().magnetic_cells {
        near(cell.parameter_work_j, 0.125 * (k1 - k0));
    }
}

#[test]
fn changed_sources_missing_regions_and_late_cancellation_publish_no_native_history() {
    let card = source::affine_magnetic_card(
        1.0 / 18.0,
        0.0,
        1.0,
        0.001,
        ElectricalPropertyBasis::Resistivity,
    );
    let problem = PlanarNonlinearMagnetostatics::new(
        common::mesh(0.3),
        materials(&card, 300.0, 8),
        source::limits(),
        &mut || false,
    )
    .unwrap();
    let original = initial(&problem);
    let at = materials(&card, 500.0, 8);
    let mut checkpoints = 0;
    original
        .step_with_materials_candidate(
            0.2,
            at.clone(),
            vec![],
            conductors(&at),
            &[9.0],
            source::limits(),
            controls(),
            &mut || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    let cancelled = original.step_with_materials_candidate(
        0.2,
        at.clone(),
        vec![],
        conductors(&at),
        &[9.0],
        source::limits(),
        controls(),
        &mut || {
            calls += 1;
            calls >= checkpoints
        },
    );
    assert!(matches!(cancelled, Err(NonlinearWindingError::Cancelled)));
    let replacement = source::affine_magnetic_card(
        1.0 / 18.0,
        0.0,
        1.0,
        0.002,
        ElectricalPropertyBasis::Resistivity,
    );
    for rejected in [materials(&replacement, 500.0, 8), at[..7].to_vec()] {
        assert!(
            original
                .step_with_materials_candidate(
                    0.2,
                    rejected,
                    vec![],
                    conductors(&at),
                    &[9.0],
                    source::limits(),
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    let hb = common::source(false, 18.0);
    let nonlinear = common::problem(&hb, 0.3);
    assert!(
        original
            .step_with_materials_candidate(
                0.2,
                nonlinear.materials().to_vec(),
                vec![],
                conductors(&at),
                &[9.0],
                source::limits(),
                controls(),
                &mut || false,
            )
            .is_err()
    );
    near(original.time_s(), 0.0);
    near(original.currents_a()[0], 1.0);
    near(original.flux_linkages_wb()[0], 0.6);
    near(original.field().potential_wb_m()[6], 0.5);
}
