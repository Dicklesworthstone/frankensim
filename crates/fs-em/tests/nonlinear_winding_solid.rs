//! Causal solid currents, native terminals and exact manufactured energy laws.
//! These numeric materials are analytical fixtures, not measured material data.

#[path = "nonlinear_winding_circuit/common.rs"]
mod common;
#[path = "nonlinear_winding_solid/slabs.rs"]
mod slabs;

use common::{controls, near};
use fs_em::{
    planar::{
        eddy::{EddyConductor, nonlinear::PlanarNonlinearEddy},
        nonlinear::PlanarNonlinearMagnetostatics,
    },
    winding::{
        nonlinear::{
            NonlinearWindingCircuit, NonlinearWindingError, NonlinearWindingInitial,
            NonlinearWindingPort, NonlinearWindings, solid::SolidWindingRegion,
        },
        ohmic::WindingConductor,
    },
};
use fs_matdb::MaterialCard;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_rep_mesh::{Metric2, TriComplex2};

const DEPTH: f64 = 0.3;
const SIGMA: f64 = 18.0;

// Two unit conducting squares have an intervening insulating square. Every
// grid vertex lies on the external trace; only the three center hats are free.
// This gives an exact nonlinear diffusion-mode oracle without joining outgoing
// and return conductors across a falsely insulated positive-area contact.
fn mesh() -> TriComplex2 {
    let mut vertices: Vec<_> = (0..=3)
        .flat_map(|x| [[f64::from(x), 0.0], [f64::from(x), 1.0]])
        .collect();
    vertices.extend((0..3).map(|x| [f64::from(x) + 0.5, 0.5]));
    let mut faces = Vec::new();
    for square in 0..3_u32 {
        let a = 2 * square;
        let center = 8 + square;
        faces.extend([
            [a, a + 2, center],
            [a + 2, a + 3, center],
            [a + 3, a + 1, center],
            [a + 1, a, center],
        ]);
    }
    TriComplex2::from_indexed_triangles(
        "solid-current-return-with-air-separation",
        vertices,
        faces,
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap()
}
fn boundary(problem: &PlanarNonlinearMagnetostatics, bias: f64) -> Vec<(usize, f64)> {
    (0..8)
        .map(|v| (v, bias * problem.mesh().vertices()[v][0]))
        .collect()
}
fn profiles() -> Vec<Vec<f64>> {
    vec![vec![
        6.0, 6.0, 6.0, 6.0, 0.0, 0.0, 0.0, 0.0, -6.0, -6.0, -6.0, -6.0,
    ]]
}
fn bindings() -> Vec<SolidWindingRegion> {
    vec![
        SolidWindingRegion {
            region: 10,
            winding: 0,
            signed_turns: 6.0,
        },
        SolidWindingRegion {
            region: 11,
            winding: 0,
            signed_turns: -6.0,
        },
    ]
}
fn sources(card: &MaterialCard, hot_triangle: Option<usize>) -> Vec<EddyConductor> {
    (0..4)
        .chain(8..12)
        .map(|triangle| EddyConductor {
            triangle,
            region: if triangle < 4 { 10 } else { 11 },
            conductor: common::material::conductor(
                card,
                if hot_triangle == Some(triangle) {
                    500.0
                } else {
                    300.0
                },
                ElectricalPropertyBasis::Resistivity,
            ),
        })
        .collect()
}
fn initial<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    core: Vec<EddyConductor>,
    air: &[usize],
    boundary: &[(usize, f64)],
) -> PlanarNonlinearEddy<'m> {
    let zero = vec![0.0; problem.mesh().faces().len()];
    let initial = problem
        .solve(
            &zero,
            boundary,
            &vec![0.0; problem.mesh().vertices().len()],
            controls().field.nonlinear,
            &mut || false,
        )
        .unwrap();
    PlanarNonlinearEddy::new_with_impressed_current(
        &initial,
        core,
        air,
        boundary,
        &zero,
        0.0,
        controls().field,
        &mut || false,
    )
    .unwrap()
}
fn port() -> NonlinearWindingPort {
    NonlinearWindingPort::VoltageInput {
        initial_voltage_v: 0.0,
        residual_scale_ohm: 1.0,
        current_tolerance_a: 1e-8,
    }
}
fn driven<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
    bias: f64,
    end: f64,
) -> NonlinearWindingCircuit<'m> {
    NonlinearWindingCircuit::new_with_solid_regions_and_ports(
        initial(
            problem,
            sources(card, None),
            &[4, 5, 6, 7],
            &boundary(problem, bias),
        ),
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 777.0 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles(),
            end_inductance_h: vec![end],
        },
        vec![port()],
        vec![],
        bindings(),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap()
}
fn mode(old: f64, current: f64, dt: f64, saturating: bool) -> f64 {
    let mut low = -2.0;
    let mut high = 2.0;
    // Integral N=1/3 and N²=1/6 on each unit square. Its projected
    // conductivity is sigma*depth/18, while H(2a) is the actual magnetic law.
    for _ in 0..80 {
        let a = f64::midpoint(low, high);
        if common::h(2.0 * a, saturating) + SIGMA * (a - old) / (36.0 * dt) > current {
            high = a;
        } else {
            low = a;
        }
    }
    f64::midpoint(low, high)
}

#[test]
fn nonlinear_pulses_retain_solid_diffusion_and_separate_physical_from_be_loss() {
    let end = 0.2;
    let resistance = 72.0 * DEPTH / SIGMA;
    for saturating in [false, true] {
        let card = common::source(saturating, SIGMA);
        let problem = common::material::problem(&card, mesh());
        for sign in [-1.0, 1.0] {
            let original = driven(&problem, &card, 0.0, end);
            let mut accepted = original.clone();
            let mut old_a = 0.0;
            let mut old_i = 0.0;
            for (target, current) in [
                (0.05, 3.0 * sign),
                (0.1, 3.0 * sign),
                (0.2, 0.0),
                (0.4, -2.0 * sign),
            ] {
                let dt = target - accepted.time_s();
                let a = mode(old_a, current, dt, saturating);
                let da = a - old_a;
                let dc_heat = resistance * current * current * dt;
                let eddy_heat = SIGMA * DEPTH * da * da / (9.0 * dt);
                let before = 2.0 * DEPTH * common::density_energy(2.0 * old_a, saturating);
                let after = 2.0 * DEPTH * common::density_energy(2.0 * a, saturating);
                let numerical = 4.0 * DEPTH * common::h(2.0 * a, saturating) * da
                    - (after - before)
                    + 0.5 * end * (current - old_i).powi(2);
                let induction = (4.0 * DEPTH * da + end * (current - old_i)) / dt;
                let voltage = induction + resistance * current;
                let candidate = accepted
                    .step_to_candidate(
                        target,
                        sources(&card, None),
                        vec![],
                        &[current],
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                let step = candidate.step();
                near(candidate.circuit().currents_a()[0], current);
                near(candidate.circuit().field().potential_wb_m()[8], a);
                near(candidate.circuit().field().potential_wb_m()[10], -a);
                near(step.inductive_voltage_v[0], induction);
                near(step.terminal_voltage_v[0], voltage);
                near(step.field.uniform_drive_work_j, dc_heat);
                near(step.field.heat_j, dc_heat + eddy_heat);
                near(step.heat_j, dc_heat + eddy_heat);
                near(step.stored_before_j, before + 0.5 * end * old_i * old_i);
                near(step.stored_after_j, after + 0.5 * end * current * current);
                near(step.numerical_dissipation_j, numerical);
                near(step.supplied_work_j, current * voltage * dt);
                assert!(step.cells.is_empty() && step.graph_resistors.is_empty());
                assert_eq!(step.field.cells.len(), 8);
                assert!(step.energy_residual_j.abs() < controls().energy_tolerance_j);
                for (region, net) in &step.field.region_currents_a {
                    near(
                        *net,
                        if *region == 10 {
                            6.0 * current
                        } else {
                            -6.0 * current
                        },
                    );
                }
                for cell in &step.field.cells {
                    let polarity = if cell.source.region == 10 { 1.0 } else { -1.0 };
                    for (local, &vertex) in problem.mesh().faces()[cell.source.triangle]
                        .iter()
                        .enumerate()
                    {
                        let change = if vertex >= 8 {
                            -2.0 * da / (3.0 * dt)
                        } else {
                            da / (3.0 * dt)
                        };
                        let electric = polarity * (6.0 * current / SIGMA + change);
                        near(cell.electric_field_v_m[local], electric);
                        near(cell.current_density_a_m2[local], SIGMA * electric);
                    }
                }
                if current == 0.0 {
                    assert!(step.field.heat_j > 0.0 && a.abs() > 0.0);
                    near(step.supplied_work_j, 0.0);
                }
                accepted = candidate.into_parts().0;
                old_a = a;
                old_i = current;
            }
            near(original.time_s(), 0.0);
            assert!(original.field().potential_wb_m().iter().all(|&a| a == 0.0));
        }
    }
}

#[test]
fn native_voltage_drive_closes_mixed_solid_stranded_and_graph_resistor_work() {
    let card = common::source(false, SIGMA);
    let problem = common::material::problem(&card, mesh());
    let copper = common::material::conductor(&card, 300.0, ElectricalPropertyBasis::Resistivity);
    let stranded: Vec<_> = [4, 6]
        .into_iter()
        .map(|triangle| WindingConductor {
            triangle,
            conductor: copper.clone(),
        })
        .collect();
    let mut profile = profiles();
    profile[0][4] = 24.0;
    profile[0][6] = -24.0;
    let end = 0.1;
    let graph_r = 0.7;
    let solid_r = 72.0 * DEPTH / SIGMA;
    let stranded_r = 288.0 * DEPTH / SIGMA;
    let mut accepted = NonlinearWindingCircuit::new_with_solid_regions(
        initial(
            &problem,
            sources(&card, None),
            &[4, 5, 6, 7],
            &boundary(&problem, 0.0),
        ),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 777.0 }),
                (2, 0, Branch::Resistor { ohms: graph_r }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profile,
            end_inductance_h: vec![end],
        },
        stranded.clone(),
        bindings(),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap();
    near(accepted.resistance_ohm()[0], solid_r + stranded_r);
    let mut old_a = 0.0;
    let mut old_i = 0.0;
    for (target, voltage) in [(0.1, 9.0), (0.3, -3.0)] {
        let dt = target - accepted.time_s();
        let q = SIGMA / (36.0 * dt);
        let coefficient = 4.0 * DEPTH / (dt * (2.0 + q));
        let current = (voltage + 2.0 * coefficient * old_a + end * old_i / dt)
            / (solid_r + stranded_r + graph_r + coefficient + end / dt);
        let a = (current + q * old_a) / (2.0 + q);
        let candidate = accepted
            .step_to_candidate(
                target,
                sources(&card, None),
                stranded.clone(),
                &[voltage],
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = candidate.step();
        near(candidate.circuit().currents_a()[0], current);
        near(candidate.circuit().field().potential_wb_m()[8], a);
        near(step.terminal_voltage_v[0], voltage - graph_r * current);
        near(
            step.cells.iter().map(|cell| cell.heat_j).sum(),
            stranded_r * current * current * dt,
        );
        near(
            step.field.uniform_drive_work_j,
            solid_r * current * current * dt,
        );
        near(
            step.graph_resistors[0].heat_j,
            graph_r * current * current * dt,
        );
        near(
            step.heat_j,
            (solid_r + stranded_r + graph_r) * current * current * dt
                + SIGMA * DEPTH * (a - old_a).powi(2) / (9.0 * dt),
        );
        near(step.supplied_work_j, voltage * current * dt);
        assert!(step.energy_residual_j.abs() < controls().energy_tolerance_j);
        accepted = candidate.into_parts().0;
        old_a = a;
        old_i = current;
    }
}

#[test]
fn conductivity_reweighting_of_stationary_bias_changes_linkage_but_never_creates_emf() {
    let card = common::source(false, SIGMA);
    let problem = common::material::problem(&card, mesh());
    let bias = 0.1;
    // Fixed affine A=bias*x carries a stationary magnetic field. Conductivity
    // is passive here and cannot supply work when its spatial weights change.
    let original = driven(&problem, &card, bias, 0.0);
    let old_potential = original.field().potential_wb_m().to_vec();
    let refreshed = original
        .step_to_candidate(
            0.2,
            sources(&card, Some(1)),
            vec![],
            &[0.0],
            controls(),
            &mut || false,
        )
        .unwrap();
    let next = refreshed.circuit();
    let step = refreshed.step();
    for (&new, &old) in next.field().potential_wb_m().iter().zip(&old_potential) {
        near(new, old);
    }
    near(next.profiles_per_m2()[0][0], 48.0 / 7.0);
    near(next.profiles_per_m2()[0][1], 24.0 / 7.0);
    near(next.resistance_ohm()[0], 540.0 * DEPTH / (7.0 * SIGMA));
    near(
        step.flux_linkages_wb[1][0] - step.flux_linkages_wb[0][0],
        -2.0 * DEPTH * bias / 7.0,
    );
    assert!((step.flux_linkages_wb[1][0] - step.flux_linkages_wb[0][0]).abs() > 1e-3);
    near(step.inductive_voltage_v[0], 0.0);
    near(step.terminal_voltage_v[0], 0.0);
    near(step.heat_j, 0.0);
    near(step.supplied_work_j, 0.0);
    near(step.numerical_dissipation_j, 0.0);
    near(step.parameter_work_j, 0.0);
    near(step.stored_after_j, step.stored_before_j);
    assert_eq!(original.profiles_per_m2(), profiles());
    assert_eq!(original.field().potential_wb_m(), old_potential);
}

#[test]
fn solid_source_refusals_and_cancelled_trials_leave_all_native_history_unchanged() {
    let card = common::source(false, SIGMA);
    let problem = common::material::problem(&card, mesh());
    let original = driven(&problem, &card, 0.0, 0.0);
    let mut shifted = sources(&card, None);
    shifted[0].region = 99;
    assert!(
        original
            .step_to_candidate(0.1, shifted, vec![], &[1.0], controls(), &mut || false)
            .is_err()
    );
    let mut limited = controls();
    limited.maximum_field_trials = 1;
    assert!(matches!(
        original.step_to_candidate(
            0.1,
            sources(&card, None),
            vec![],
            &[1.0],
            limited,
            &mut || false,
        ),
        Err(NonlinearWindingError::TrialBudget)
    ));
    assert!(matches!(
        original.step_to_candidate(
            0.1,
            sources(&card, None),
            vec![],
            &[1.0],
            controls(),
            &mut || true,
        ),
        Err(NonlinearWindingError::Cancelled)
    ));
    near(original.time_s(), 0.0);
    near(original.currents_a()[0], 0.0);
    near(original.flux_linkages_wb()[0], 0.0);
    assert!(original.field().potential_wb_m().iter().all(|&a| a == 0.0));
}
