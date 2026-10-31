//! A spatial continuum oracle for the first causal diffusion interval.
use super::*;

const HEIGHT: f64 = 0.25;
const DEPTH: f64 = 0.4;
const SIGMA: f64 = 2.0;
const DT: f64 = 0.125;
const STARTS: [f64; 4] = [-3.0, -1.0, 1.0, 3.0];

fn bar(x: f64) -> Option<usize> {
    STARTS
        .iter()
        .position(|&start| x > start && x < start + 1.0)
}
fn slab_mesh(h: f64) -> TriComplex2 {
    let nx = (9.0 / h) as usize;
    let ny = (HEIGHT / h) as usize;
    let vertices = (0..=nx)
        .flat_map(|x| (0..=ny).map(move |y| [-4.0 + x as f64 * h, y as f64 * h]))
        .collect();
    let mut faces = Vec::new();
    for x in 0..nx {
        for y in 0..ny {
            let a = (x * (ny + 1) + y) as u32;
            let b = a + 1;
            let c = a + (ny + 1) as u32;
            faces.extend([[a, c, c + 1], [a, c + 1, b]]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "causal-four-solid-slab-current-constraints",
        vertices,
        faces,
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap()
}
fn model<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
) -> NonlinearWindingCircuit<'m> {
    let mut profiles = vec![vec![0.0; problem.mesh().faces().len()]; 2];
    let mut core = Vec::new();
    let mut air = Vec::new();
    let copper = common::material::conductor(card, 300.0, ElectricalPropertyBasis::Resistivity);
    for (triangle, face) in problem.mesh().faces().iter().enumerate() {
        let x = face
            .iter()
            .map(|&v| problem.mesh().vertices()[v as usize][0] / 3.0)
            .sum();
        if let Some(section) = bar(x) {
            let winding = if section == 0 || section == 3 { 0 } else { 1 };
            profiles[winding][triangle] = if section < 2 {
                HEIGHT.recip()
            } else {
                -HEIGHT.recip()
            };
            core.push(EddyConductor {
                triangle,
                region: 10 + section as u32,
                conductor: copper.clone(),
            });
        } else {
            air.push(triangle);
        }
    }
    let boundary: Vec<_> = problem
        .mesh()
        .vertices()
        .iter()
        .enumerate()
        .filter(|(_, point)| point[0] == -4.0)
        .map(|(v, _)| (v, 0.0))
        .collect();
    let bindings = (0..4)
        .map(|section| SolidWindingRegion {
            region: 10 + section,
            winding: if section == 0 || section == 3 { 0 } else { 1 },
            signed_turns: if section < 2 { 1.0 } else { -1.0 },
        })
        .collect();
    NonlinearWindingCircuit::new_with_solid_regions_and_ports(
        initial(problem, core, &air, &boundary),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 777.0 }),
                (2, 0, Branch::CurrentSource { port: 1 }),
                (2, 0, Branch::Inductor { henries: 999.0 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1, 3],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.01, 0.02],
        },
        vec![port(), port()],
        vec![],
        bindings,
        NonlinearWindingInitial {
            currents_a: vec![0.0; 2],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0; 2],
        },
        common::config(),
        &mut || false,
    )
    .unwrap()
}

// At the first backward-Euler interval from A=0, J''=(mu*sigma/dt)J.
// Interface derivatives are gamma² times the cumulative axial current/HEIGHT.
// This continuous slab solution does not use a native field matrix, recovered
// A, circuit solve, fitted inductance or frequency-domain approximation.
fn analytic_current(section: usize, x: f64, currents: [f64; 2]) -> f64 {
    let regional = [currents[0], currents[1], -currents[1], -currents[0]];
    let before: f64 = regional[..section].iter().sum();
    let after = before + regional[section];
    let gamma = (SIGMA / DT).sqrt();
    let s = gamma * (x - STARTS[section]);
    gamma / HEIGHT * ((after - before * gamma.cosh()) * s.cosh() / gamma.sinh() + before * s.sinh())
}
fn analytic_power(currents: [f64; 2]) -> f64 {
    let mut integral = 0.0;
    let steps = 256;
    for section in 0..4 {
        for sample in 0..=steps {
            let weight = if sample == 0 || sample == steps {
                1.0
            } else if sample % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let j = analytic_current(
                section,
                STARTS[section] + sample as f64 / steps as f64,
                currents,
            );
            integral += weight * j * j / (3.0 * steps as f64);
        }
    }
    HEIGHT * DEPTH * integral / SIGMA
}

#[test]
fn two_winding_skin_and_proximity_profiles_converge_to_causal_slab_diffusion() {
    let card = common::source(false, SIGMA);
    let currents = [0.018, 0.013];
    let expected_heat = DT * analytic_power(currents);
    let mut errors = Vec::new();
    for h in [0.25, 0.125] {
        let problem = common::material::problem(&card, slab_mesh(h));
        let original = model(&problem, &card);
        let candidate = original
            .step_to_candidate(
                DT,
                original.field().conductors().to_vec(),
                vec![],
                &currents,
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = candidate.step();
        let mut square_error = 0.0;
        let mut square_reference = 0.0;
        let mut redistribution = 0.0_f64;
        let mut local_cross_term = 0.0_f64;
        for cell in &step.field.cells {
            let section = (cell.source.region - 10) as usize;
            let winding = if section == 0 || section == 3 { 0 } else { 1 };
            let dc = currents[winding] / HEIGHT * if section < 2 { 1.0 } else { -1.0 };
            for (local, &v) in problem.mesh().faces()[cell.source.triangle]
                .iter()
                .enumerate()
            {
                let expected =
                    analytic_current(section, problem.mesh().vertices()[v as usize][0], currents);
                let actual = cell.current_density_a_m2[local];
                square_error += cell.volume_m3 * (actual - expected).powi(2) / 3.0;
                square_reference += cell.volume_m3 * expected * expected / 3.0;
                redistribution = redistribution.max((actual - dc).abs());
                near(actual, SIGMA * cell.electric_field_v_m[local]);
            }
            let eddy = cell.current_density_a_m2.map(|j| j - dc);
            let sum: f64 = eddy.iter().sum();
            let separated = cell.volume_m3 / (12.0 * SIGMA)
                * (sum * sum + eddy.iter().map(|j| j * j).sum::<f64>())
                + cell.volume_m3 * dc * dc / SIGMA;
            local_cross_term = local_cross_term.max((cell.power_w - separated).abs());
        }
        errors.push((square_error / square_reference).sqrt());
        assert!(redistribution > 0.02);
        assert!(local_cross_term > 1e-8);
        assert!((step.field.heat_j - expected_heat).abs() / expected_heat < 0.15);
        near(
            step.field.uniform_drive_work_j,
            DT * 2.0 * DEPTH / (SIGMA * HEIGHT) * currents.iter().map(|i| i * i).sum::<f64>(),
        );
        near(step.heat_j, step.field.heat_j);
        assert!(step.cells.is_empty() && step.graph_resistors.is_empty());
        assert!(step.energy_residual_j.abs() < controls().energy_tolerance_j);
        for &(region, current) in &step.field.region_currents_a {
            let section = (region - 10) as usize;
            near(
                current,
                [currents[0], currents[1], -currents[1], -currents[0]][section],
            );
        }
        if h == 0.125 {
            // The inner terminal current is held fixed. Reversing only its
            // neighbor changes the inner local current profile through the
            // same causal A history, demonstrating actual proximity coupling.
            let changed = original
                .step_to_candidate(
                    DT,
                    original.field().conductors().to_vec(),
                    vec![],
                    &[-currents[0], currents[1]],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            let largest = step
                .field
                .cells
                .iter()
                .zip(&changed.step().field.cells)
                .filter(|(cell, _)| cell.source.region == 11)
                .flat_map(|(a, b)| {
                    a.current_density_a_m2
                        .into_iter()
                        .zip(b.current_density_a_m2)
                        .map(|(a, b)| (a - b).abs())
                })
                .fold(0.0_f64, f64::max);
            assert!(largest > 0.02);
            near(changed.circuit().currents_a()[1], currents[1]);
        }
    }
    assert!(errors[1] < 0.8 * errors[0], "spatial errors: {errors:?}");
    assert!(errors[1] < 0.08, "fine spatial error: {}", errors[1]);
}
