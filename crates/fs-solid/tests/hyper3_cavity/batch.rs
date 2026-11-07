//! Shared full-mesh admission must retain every cavity's original local chart.

use super::*;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_solid::hyper3::pressure::cavity::{AdiabaticCavity, TetCavityGasLoad};

fn same_bits(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(actual.to_bits(), expected.to_bits(), "coordinate {index}");
    }
}

fn same_load(
    actual: &TetCavityGasLoad,
    expected: &TetCavityGasLoad,
    direction: &[f64],
    cx: &Cx<'_>,
) {
    assert_eq!(actual.gas, expected.gas);
    assert_eq!(actual.volume_m3.to_bits(), expected.volume_m3.to_bits());
    assert_eq!(
        actual.net_pressure_pa.to_bits(),
        expected.net_pressure_pa.to_bits()
    );
    assert_eq!(
        actual.potential_change_j.to_bits(),
        expected.potential_change_j.to_bits()
    );
    same_bits(&actual.force_n, &expected.force_n);
    same_bits(
        actual.force_tangent.acoustic_rank_one_vector(),
        expected.force_tangent.acoustic_rank_one_vector(),
    );
    same_bits(
        &actual.force_tangent.apply(direction, cx).unwrap(),
        &expected.force_tangent.apply(direction, cx).unwrap(),
    );
}

#[test]
fn g3_batch_cavities_preserve_tiny_local_motion_across_remote_reference_charts() {
    with_cx(|_, cx| {
        let (local, local_tets) = shell();
        let far = 2.0_f64.powi(36);
        let offset = [far, -far, far];
        let mut nodes = local.clone();
        nodes.extend(
            local
                .iter()
                .map(|x| std::array::from_fn(|axis| x[axis] + offset[axis])),
        );
        let mut tets = local_tets.clone();
        tets.extend(
            local_tets
                .iter()
                .map(|tet| tet.map(|node| node + local.len())),
        );

        // These separately prepared surfaces are equal reference charts, not
        // aliases of one Arc. Each source owns its own selected inner boundary.
        let surfaces = [surface(&nodes, &tets, cx), surface(&nodes, &tets, cx)];
        assert!(!Arc::ptr_eq(&surfaces[0], &surfaces[1]));
        let cavities: [_; 2] = std::array::from_fn(|i| {
            let first = i * local.len() + 4;
            let faces: Vec<_> = surfaces[i]
                .faces()
                .iter()
                .enumerate()
                .filter_map(|(face, nodes)| {
                    nodes
                        .iter()
                        .all(|node| (first..first + 4).contains(node))
                        .then_some(face)
                })
                .collect();
            ClosedTetCavity::new(Arc::clone(&surfaces[i]), &faces, cx).unwrap()
        });
        let gases = [
            super::gas_fixture::gas(cavities[0].reference_volume_m3(), 1e5, 0.001, 0.0),
            super::gas_fixture::gas(cavities[1].reference_volume_m3(), 2e5, 0.002, 1e24),
        ];
        let controls = GasEntropySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 600.0,
            entropy_tolerance_j_kg_k: 1e-10,
            maximum_iterations: 128,
        };
        let sources: [_; 2] = std::array::from_fn(|i| AdiabaticCavity {
            cavity: &cavities[i],
            reference: &gases[i],
            exterior_pressure_pa: 5e4,
            controls,
        });
        let common = [2.0_f64.powi(20), -2.0_f64.powi(19), 2.0_f64.powi(18)];
        let tiny = 2.0_f64.powi(-23);
        let mut displacement = Vec::with_capacity(3 * nodes.len());
        for chamber in 0..2 {
            let stretch_increment = if chamber == 0 { -0.03 } else { tiny };
            for x in &local {
                for axis in 0..3 {
                    displacement.push(common[axis] + stretch_increment * x[axis]);
                }
            }
        }
        let expected = sources
            .each_ref()
            .map(|source| source.evaluate(&displacement, cx).unwrap());

        // The small displacement is represented; it vanishes only if someone
        // adds reference and displacement globally BEFORE taking local offsets.
        let global_inner = std::array::from_fn(|i| {
            std::array::from_fn(|axis| {
                let node = local.len() + 4 + i;
                nodes[node][axis] + displacement[3 * node + axis]
            })
        });
        let rounded_global_volume = tetra_volume_gradient(global_inner).0;
        assert_eq!(
            rounded_global_volume.to_bits(),
            cavities[1].reference_volume_m3().to_bits()
        );
        assert!(expected[1].volume_m3 > cavities[1].reference_volume_m3());
        assert_ne!(
            expected[1].volume_m3.to_bits(),
            rounded_global_volume.to_bits()
        );
        assert_ne!(expected[1].potential_change_j, 0.0);

        let prescribed: Vec<_> = displacement.iter().copied().enumerate().collect();
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2e6,
                lambda: 3e6,
            },
            4.0,
        )
        .unwrap();
        let problem = HyperTetProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            material: HyperelasticLawField::Uniform(&material),
            prescribed_m: &prescribed,
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings {
                load_steps: 1,
                ..HyperTetSettings::default()
            },
        };
        let direction: Vec<_> = (0..displacement.len())
            .map(|i| ((i * 7 % 19) as f64 - 9.0) / 16.0)
            .collect();
        let mut original_solid = None;
        for order in [[0, 1], [1, 0]] {
            let ordered = order.map(|i| sources[i]);
            let solution = problem
                .solve_with_adiabatic_cavities(&ordered, None, None, cx)
                .unwrap();
            let actual: Vec<_> = solution
                .solid
                .displacement_m
                .iter()
                .flatten()
                .copied()
                .collect();
            same_bits(&actual, &displacement);
            for (load, index) in solution.cavities.iter().zip(order) {
                same_load(load, &expected[index], &direction, cx);
            }
            if let Some(original) = &original_solid {
                assert_eq!(&solution.solid, original);
            } else {
                original_solid = Some(solution.solid);
            }
        }
    });
}
