//! Oriented compressive plastic flow, reversal and confinement through native FE.
//!
//! Coefficients are synthetic. A scalar uniaxial return supplies independent
//! stress and flow oracles, including an oblique material frame with plastic
//! shear. Lateral unknowns remain free except for rigid-motion constraints.
//! Hill history stays in material coordinates; the observations below are
//! spatial quantities and use the original orientation owner's conversion.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::IsotropicElastic;
use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::{Hill48, Hill48Coefficients, HillState};
use fs_material::plastic::reduction::ReductionConfig;
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem, SmallStrainPlaneSolution, SmallStrainTetProblem,
    SmallStrainTetSolution, TetAssemblyBudget,
};

const E: f64 = 70e9;
const NU: f64 = 0.3;
const YIELD: f64 = 100e6;
const HARDENING: f64 = 1e9;
const COEFFICIENTS: Hill48Coefficients = Hill48Coefficients {
    f: 0.2,
    g: 1.2,
    h: 0.6,
    l: 1.8,
    m: 1.6,
    n: 1.1,
};
const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [2.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [2.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    [2.0, 1.0, 1.0],
];
const TETS: [[usize; 4]; 6] = [
    [0, 1, 3, 7],
    [0, 1, 7, 5],
    [0, 2, 7, 3],
    [0, 2, 6, 7],
    [0, 4, 5, 7],
    [0, 4, 7, 6],
];

#[derive(Clone, Copy)]
struct Direction {
    rotation: [[f64; 3]; 3],
    squared_factor: f64,
    /// The spatial associated flow normal is sign(sigma) times this / a.
    flow_numerator: [f64; 6],
}

const C: f64 = std::f64::consts::FRAC_1_SQRT_2;
const DIRECTIONS: [Direction; 4] = [
    Direction {
        rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        squared_factor: 1.8, // G + H; material x follows spatial x.
        flow_numerator: [1.8, -0.6, -1.2, 0.0, 0.0, 0.0],
    },
    Direction {
        rotation: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        squared_factor: 0.8, // F + H; material y follows spatial x.
        flow_numerator: [0.8, -0.6, -0.2, 0.0, 0.0, 0.0],
    },
    Direction {
        rotation: [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        squared_factor: 1.4, // F + G; material z follows spatial x.
        flow_numerator: [1.4, -1.2, -0.2, 0.0, 0.0, 0.0],
    },
    Direction {
        rotation: [[C, -C, 0.0], [C, C, 0.0], [0.0, 0.0, 1.0]],
        // At 45 degrees, local sigma11=sigma22=sigma/2 and sigma12=-sigma/2.
        // Thus a²=(F+G+2N)/4. Rotating the associated flow back gives the
        // spatial numerators [a²,(F+G-2N)/4,-(F+G)/2,(G-F)/4,0,0].
        squared_factor: 0.9,
        flow_numerator: [0.9, -0.2, -0.7, 0.25, 0.0, 0.0],
    },
];

fn material(direction: Direction) -> OrientedHill48 {
    let local = Hill48::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        YIELD,
        HARDENING,
        COEFFICIENTS,
    )
    .unwrap();
    OrientedHill48::new(local, direction.rotation).unwrap()
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "actual {actual:.15e}, expected {expected:.15e}, tolerance {tolerance:e}"
    );
}

#[derive(Default)]
struct AxialOracle {
    plastic: [f64; 6],
    alpha: f64,
}

impl AxialOracle {
    fn advance(&mut self, strain: f64, direction: Direction) -> f64 {
        // Uniaxial stress gives q = a |sigma| and dp_x = sign(sigma) a dgamma.
        // Hooke's law and q = Y + H alpha yield E a² + H independently of
        // the tensor return. The same scalar formula follows load reversal.
        let a = direction.squared_factor.sqrt();
        let trial = E * (strain - self.plastic[0]);
        let gamma = (a * trial.abs() - YIELD - HARDENING * self.alpha).max(0.0)
            / (E * direction.squared_factor + HARDENING);
        for (plastic, numerator) in self.plastic.iter_mut().zip(direction.flow_numerator) {
            *plastic += trial.signum() * gamma * numerator / a;
        }
        self.alpha += gamma;
        E * (strain - self.plastic[0])
    }

    fn strain(&self, axial: f64, sigma: f64) -> [f64; 6] {
        let mut strain = self.plastic;
        strain[0] = axial;
        strain[1] -= NU * sigma / E;
        strain[2] -= NU * sigma / E;
        strain
    }

    fn check(&self, law: &OrientedHill48, state: &HillState) {
        let spatial = law.plastic_strain_spatial(state);
        for (&actual, &expected) in spatial.iter().zip(&self.plastic) {
            close(actual, expected, 2e-9);
        }
        close(state.alpha, self.alpha, 2e-9);
        close(spatial[..3].iter().sum(), 0.0, 1e-11);
    }
}

fn tet_supports() -> Vec<(usize, f64)> {
    let mut supports = Vec::new();
    for (node, x) in NODES.iter().enumerate() {
        supports.push((3 * node, x[0]));
        if *x == [0.0, 0.0, 0.0] {
            supports.extend([(3 * node + 1, 0.0), (3 * node + 2, 0.0)]);
        }
        if *x == [0.0, 1.0, 0.0] {
            supports.push((3 * node + 2, 0.0));
        }
    }
    // All nodes lie on end faces. Only rigid modes are pinned laterally;
    // plastic xy shear can produce u_y = 2 eps_xy x + eps_yy y.
    assert!(supports.len() < 3 * NODES.len());
    supports
}

fn tet_problem<'a>(
    law: &'a OrientedHill48,
    supports: &'a [(usize, f64)],
    path: &'a [f64],
) -> SmallStrainTetProblem<'a, OrientedHill48> {
    SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(law),
        prescribed_m: supports,
        nodal_forces_n: &[],
        load_path: path,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    }
}

#[test]
fn rotated_tet_compression_and_reversal_resolve_free_transverse_and_shear_motion() {
    let supports = tet_supports();
    for direction in DIRECTIONS {
        let law = material(direction);
        let mut accepted: Option<SmallStrainTetSolution<HillState>> = None;
        let mut oracle = AxialOracle::default();
        for strain in [-0.0005, -0.0015, -0.004, -0.001, 0.004] {
            let path = [strain];
            let problem = tet_problem(&law, &supports, &path);
            let next = with_cx(|cx| match &accepted {
                None => problem.solve(cx),
                Some(old) => {
                    problem.solve_from(old.displacement_m.as_flattened(), &old.element_states, cx)
                }
            })
            .unwrap();
            let sigma = oracle.advance(strain, direction);
            let expected_strain = oracle.strain(strain, sigma);
            for element in 0..TETS.len() {
                oracle.check(&law, &next.element_states[element]);
                for component in 0..6 {
                    close(
                        next.element_strain[element][component],
                        expected_strain[component],
                        2e-9,
                    );
                    close(
                        next.element_stress[element][component],
                        if component == 0 { sigma } else { 0.0 },
                        5.0,
                    );
                }
            }
            for (x, u) in NODES.iter().zip(&next.displacement_m) {
                close(u[0], strain * x[0], 2e-9);
                close(
                    u[1],
                    expected_strain[1] * x[1] + 2.0 * expected_strain[3] * x[0],
                    6e-9,
                );
                close(u[2], expected_strain[2] * x[2], 2e-9);
            }
            let reaction: f64 = supports
                .iter()
                .zip(&next.increments[0].reactions_n)
                .filter(|((dof, _), _)| dof % 3 == 0 && NODES[dof / 3][0] == 2.0)
                .map(|(_, force)| force)
                .sum();
            close(reaction, sigma, 5.0); // Loaded area is 1 m².
            if strain == -0.0015 {
                assert_eq!(
                    next.increments[0].evolving_elements > 0,
                    E * strain.abs() * direction.squared_factor.sqrt() > YIELD
                );
            }
            accepted = Some(next);
        }
    }
}

fn plane_supports(mesh: &Mesh2, driven: bool) -> Vec<(usize, f64)> {
    let mut supports = Vec::new();
    for (node, x) in mesh.nodes.iter().enumerate() {
        if x[0] == 0.0 || (driven && x[0] == 2.0) {
            supports.push((2 * node, if driven { x[0] } else { 0.0 }));
        }
        if *x == [0.0, 0.0] {
            supports.push((2 * node + 1, 0.0));
        }
    }
    assert!(supports.len() < 2 * mesh.nodes.len());
    supports
}

fn plane_problem<'a>(
    mesh: &'a Mesh2,
    law: &'a OrientedHill48,
    supports: &'a [(usize, f64)],
    path: &'a [f64],
) -> SmallStrainPlaneProblem<'a, OrientedHill48> {
    SmallStrainPlaneProblem {
        mesh,
        laws: SmallStrainLawField::Uniform(law),
        plane: PlaneKind::Stress,
        formulation: Formulation::Standard,
        thickness_m: 0.2,
        reduction: ReductionConfig {
            max_iterations: 60,
            stress_tolerance_pa: 1e-4,
            strain_limit: 0.02,
        },
        rotation_limit_rad: 0.02,
        prescribed_m: supports,
        nodal_forces_n: &[],
        load_path: path,
        budget: PlaneAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    }
}

#[test]
fn rotated_p1_q1_compression_resolves_width_thickness_shear_and_free_release() {
    for mesh in [
        Mesh2::quads(2.0, 1.0, 3, 2),
        Mesh2::triangles(2.0, 1.0, 3, 2),
    ] {
        let supports = plane_supports(&mesh, true);
        // One out-of-plane permutation and one oblique in-plane frame.
        for direction in [DIRECTIONS[2], DIRECTIONS[3]] {
            let law = material(direction);
            let mut oracle = AxialOracle::default();
            let mut accepted: Option<SmallStrainPlaneSolution<HillState>> = None;
            for strain in [-0.0005, -0.0015, -0.004] {
                let path = [strain];
                let problem = plane_problem(&mesh, &law, &supports, &path);
                let next = with_cx(|cx| match &accepted {
                    None => problem.solve(cx),
                    Some(old) => problem.solve_from(
                        old.displacement_m.as_flattened(),
                        &old.quadrature_states,
                        cx,
                    ),
                })
                .unwrap();
                let sigma = oracle.advance(strain, direction);
                let expected_strain = oracle.strain(strain, sigma);
                for (state, stress) in next
                    .quadrature_states
                    .iter()
                    .zip(&next.quadrature_stress_pa)
                {
                    oracle.check(&law, &state.history);
                    for component in 0..6 {
                        close(state.strain[component], expected_strain[component], 2e-9);
                        close(
                            stress[component],
                            if component == 0 { sigma } else { 0.0 },
                            5.0,
                        );
                    }
                }
                let reaction: f64 = supports
                    .iter()
                    .zip(&next.increments[0].reactions_n)
                    .filter(|((dof, _), _)| dof % 2 == 0 && mesh.nodes[dof / 2][0] == 2.0)
                    .map(|(_, force)| force)
                    .sum();
                close(reaction, 0.2 * sigma, 1.0);
                accepted = Some(next);
            }
            let compressed = accepted.unwrap();
            let preserved = compressed.clone();
            let released = plane_supports(&mesh, false);
            let unloaded = with_cx(|cx| {
                plane_problem(&mesh, &law, &released, &[0.0]).solve_from(
                    compressed.displacement_m.as_flattened(),
                    &compressed.quadrature_states,
                    cx,
                )
            })
            .unwrap();
            assert_eq!(compressed, preserved);
            assert!(oracle.plastic[0] < 0.0 && oracle.plastic[1] > 0.0 && oracle.plastic[2] > 0.0);
            for (x, u) in mesh.nodes.iter().zip(&unloaded.displacement_m) {
                close(u[0], oracle.plastic[0] * x[0], 4e-9);
                close(
                    u[1],
                    oracle.plastic[1] * x[1] + 2.0 * oracle.plastic[3] * x[0],
                    6e-9,
                );
            }
            for (state, stress) in unloaded
                .quadrature_states
                .iter()
                .zip(&unloaded.quadrature_stress_pa)
            {
                oracle.check(&law, &state.history);
                for component in 0..6 {
                    close(state.strain[component], oracle.plastic[component], 2e-9);
                    close(stress[component], 0.0, 5.0);
                }
            }
        }
    }
}

#[test]
fn oriented_q1_bbar_retains_plane_strain_confinement_and_material_history() {
    let mesh = Mesh2::quads(2.0, 1.0, 3, 2);
    let law = material(DIRECTIONS[3]);
    let supports = plane_supports(&mesh, true);
    let mut standard = plane_problem(&mesh, &law, &supports, &[-0.0005, -0.002, -0.004]);
    standard.plane = PlaneKind::Strain;
    let mut bbar = standard.clone();
    bbar.formulation = Formulation::BBar;
    let compatible = with_cx(|cx| standard.solve(cx)).unwrap();
    let projected = with_cx(|cx| bbar.solve(cx)).unwrap();
    // Constant volumetric strain is preserved by B-bar. The genuine lateral
    // unknowns distinguish this from an all-prescribed patch; this is not a
    // volumetric-locking or mesh-convergence benchmark.
    for (left, right) in compatible
        .displacement_m
        .iter()
        .zip(&projected.displacement_m)
    {
        for component in 0..2 {
            close(left[component], right[component], 2e-9);
        }
    }
    for ((state, expected), (stress, expected_stress)) in projected
        .quadrature_states
        .iter()
        .zip(&compatible.quadrature_states)
        .zip(
            projected
                .quadrature_stress_pa
                .iter()
                .zip(&compatible.quadrature_stress_pa),
        )
    {
        let plastic = law.plastic_strain_spatial(&state.history);
        assert!(state.history.alpha > 0.0);
        close(state.strain[2], 0.0, 2e-11);
        close(plastic[..3].iter().sum(), 0.0, 1e-11);
        close(stress[1], 0.0, 5.0);
        close(stress[3], 0.0, 5.0);
        assert!(stress[2] < -1e6);
        assert!(plastic[2] > 0.0);
        for component in 0..6 {
            close(
                state.history.plastic_strain[component],
                expected.history.plastic_strain[component],
                2e-9,
            );
            close(stress[component], expected_stress[component], 5.0);
        }
    }
}
