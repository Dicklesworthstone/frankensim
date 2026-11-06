//! Actual native solid equilibrium with the original finite gas parcel.
use super::*;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::{HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField};

fn material() -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu: 200_000.0, lambda: 300_000.0 }, 4.0).unwrap()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]], tets: &'a [[usize; 4]], law: &'a Hyperelastic,
    supports: &'a [(usize, f64)], forces: &'a [(usize, f64)], steps: usize,
) -> HyperTetProblem<'a> {
    HyperTetProblem {
        nodes_m: nodes, tetrahedra: tets, material: HyperelasticLawField::Uniform(law),
        prescribed_m: supports, nodal_forces_n: forces, body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings { load_steps: steps, max_newton_iterations: 60,
            relative_tolerance: 1e-9, ..HyperTetSettings::default() },
    }
}

#[test]
fn sealed_cavity_native_equilibrium_keeps_source_feedback_and_reaction_balance() {
    with_cx(|_, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let law = material();
        let supports: Vec<_> = (0..12).map(|dof| (dof, 0.0)).collect();
        let forces = [(12, 1000.0)];
        for slope in [0.0, 0.001] {
            let seed = gas(0.75, slope);
            let reference = ClosedGasParcel::new(seed.model().clone(), seed.point().clone(),
                cavity.reference_volume_m3()).unwrap();
            let mut first: Option<Vec<f64>> = None;
            for steps in [1, 4] {
                let problem = problem(&nodes, &tets, &law, &supports, &forces, steps);
                let result = problem.solve_with_adiabatic_cavity(&cavity, &reference, 30_000.0, entropy(), cx).unwrap();
                let displacement: Vec<_> = result.displacement_m.iter().flatten().copied().collect();
                let endpoint = cavity.adiabatic_gas_load(&displacement, &reference, 30_000.0, entropy(), cx).unwrap();
                assert!(endpoint.volume_m3 > reference.volume_m3());
                assert!(endpoint.gas.state().state().pressure < reference.state().state().pressure);
                close(endpoint.gas.mass_kg(), reference.mass_kg(), 2e-13);
                let internal = problem.internal_force(&displacement, cx).unwrap();
                let dead = problem.external_force_n(cx).unwrap();
                let mut force_scale = 0.0_f64;
                let mut residual = 0.0_f64;
                for dof in 12..24 {
                    let external = dead[dof] + endpoint.force_n[dof];
                    force_scale = force_scale.hypot(external);
                    residual = residual.hypot(internal[dof] - external);
                }
                assert!(residual <= 5e-8 * force_scale, "{residual:e} / {force_scale:e}");
                let mut reaction_sum = [0.0; 3];
                for &(dof, reaction) in &result.reactions_n {
                    close(reaction, internal[dof] - dead[dof] - endpoint.force_n[dof], 1e-10);
                    reaction_sum[dof % 3] += reaction;
                }
                for (axis, expected) in [-1000.0, 0.0, 0.0].into_iter().enumerate() {
                    assert!((reaction_sum[axis] - expected).abs() < 0.005);
                }
                if slope == 0.0 {
                    // Independent ideal-gas isentrope, using the source's
                    // constant mixture gamma only in this reference check.
                    let initial = reference.state().state();
                    let expected = initial.pressure
                        * (reference.volume_m3() / endpoint.volume_m3).powf(initial.gamma);
                    close(endpoint.gas.state().state().pressure, expected, 2e-10);
                }
                if let Some(previous) = &first {
                    for (&a, &b) in displacement.iter().zip(previous) { close(a, b, 2e-8); }
                } else { first = Some(displacement); }
            }
            close(reference.state().state().pressure, 100_000.0, 1e-13);
            close(reference.state().state().temperature, 400.0, 1e-13);
        }
    });
}

#[test]
fn sealed_cavity_resists_prescribed_compression_through_its_native_acoustic_stiffness() {
    with_cx(|_, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let seed = gas(0.75, 0.001);
        let reference = ClosedGasParcel::new(seed.model().clone(), seed.point().clone(),
            cavity.reference_volume_m3()).unwrap();
        let mut supports: Vec<_> = nodes[..4].iter().enumerate().flat_map(|(n, x)|
            (0..3).map(move |axis| (3*n+axis, -0.04*x[axis]))).collect();
        // This inner-wall coordinate has a nonzero volume derivative, so the
        // pending motion exercises the off-graph acoustic K_fp contribution.
        supports.push((15, -0.04 * nodes[5][0]));
        let law = material();
        let problem = problem(&nodes, &tets, &law, &supports, &[], 4);
        let result = problem.solve_with_adiabatic_cavity(&cavity, &reference, 100_000.0, entropy(), cx).unwrap();
        let displacement: Vec<_> = result.displacement_m.iter().flatten().copied().collect();
        let endpoint = cavity.adiabatic_gas_load(&displacement, &reference, 100_000.0, entropy(), cx).unwrap();
        assert!(endpoint.volume_m3 < reference.volume_m3());
        assert!(endpoint.gas.state().state().pressure > reference.state().state().pressure);
        assert!(endpoint.gas.state().state().temperature > reference.state().state().temperature);
        let empty = problem.solve(cx).unwrap();
        let empty_displacement: Vec<_> = empty.displacement_m.iter().flatten().copied().collect();
        assert!(endpoint.volume_m3 > cavity.evaluate(&empty_displacement, cx).unwrap().volume_m3());
        for &(dof, expected) in &supports { close(displacement[dof], expected, 1e-13); }
        let internal = problem.internal_force(&displacement, cx).unwrap();
        let residual = (12..24).filter(|&dof| dof != 15)
            .fold(0.0_f64, |r, dof| r.hypot(internal[dof] - endpoint.force_n[dof]));
        let scale = (12..24).filter(|&dof| dof != 15)
            .fold(0.0_f64, |r, dof| r.hypot(endpoint.force_n[dof]));
        assert!(residual <= 5e-8 * scale, "{residual:e} / {scale:e}");
        let reaction = result.reactions_n.iter().find(|&&(dof, _)| dof == 15).unwrap().1;
        close(reaction, internal[15] - endpoint.force_n[15], 1e-10);
    });
}

#[test]
fn sealed_cavity_native_admission_and_cancellation_preserve_the_reference_parcel() {
    with_cx(|gate, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let seed = gas(0.75, 0.001);
        let reference = ClosedGasParcel::new(seed.model().clone(), seed.point().clone(),
            cavity.reference_volume_m3()).unwrap();
        let law = material();
        let supports: Vec<_> = (0..12).map(|dof| (dof, 0.0)).collect();
        let problem = problem(&nodes, &tets, &law, &supports, &[], 1);
        assert!(problem.solve_with_adiabatic_cavity(&cavity, &seed, 30_000.0, entropy(), cx).is_err());
        assert!(problem.solve_with_adiabatic_cavity(&cavity, &reference, -1.0, entropy(), cx).is_err());
        let mut shifted = nodes.clone();
        shifted[0][0] += 0.01;
        let other = HyperTetProblem { nodes_m: &shifted, ..problem };
        assert!(other.solve_with_adiabatic_cavity(&cavity, &reference, 30_000.0, entropy(), cx).is_err());
        gate.request();
        assert!(matches!(problem.solve_with_adiabatic_cavity(&cavity, &reference, 30_000.0, entropy(), cx),
            Err(HyperTetError::Cancelled)));
        close(reference.volume_m3(), cavity.reference_volume_m3(), 1e-15);
        close(reference.state().state().pressure, 100_000.0, 1e-13);
    });
}
