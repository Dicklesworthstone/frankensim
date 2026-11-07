//! Distinct sealed inventories coupled through one connected native solid.
use super::gas_fixture::gas;
use super::*;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_solid::hyper3::gravity::{GravitySample3, ReferenceGravityDensity, StationaryGravity3};
use fs_solid::hyper3::pressure::cavity::AdiabaticCavity;

const FIXED: [(usize, f64); 6] = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];

fn controls() -> GasEntropySolve {
    GasEntropySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 600.0,
        entropy_tolerance_j_kg_k: 1e-10,
        maximum_iterations: 128,
    }
}

/// Reflect the existing conforming hollow tetrahedron across its outer z=0
/// face. Shared vertices and the shared triangle are genuine solid connectivity;
/// only the two inner boundaries remain separate pressure surfaces.
pub(super) fn connected_shells() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let (mut nodes, mut tets) = shell();
    let reflected = [0, 1, 2, 8, 9, 10, 11, 12];
    for node in [3, 4, 5, 6, 7] {
        let [x, y, z] = nodes[node];
        nodes.push([x, y, -z]);
    }
    for index in 0..12 {
        let mut tet = tets[index].map(|node| reflected[node]);
        tet.swap(0, 1); // Reflection reverses the determinant.
        tets.push(tet);
    }
    assert_eq!((nodes.len(), tets.len()), (13, 24));
    (nodes, tets)
}

pub(super) fn cavities(nodes: &[[f64; 3]], tets: &[[usize; 4]], cx: &Cx<'_>) -> [ClosedTetCavity; 2] {
    let surface = surface(nodes, tets, cx);
    assert_eq!(
        surface.faces().len(),
        14,
        "the joining outer face must disappear"
    );
    std::array::from_fn(|chamber| {
        let first = if chamber == 0 { 4 } else { 9 };
        let faces: Vec<_> = surface
            .faces()
            .iter()
            .enumerate()
            .filter_map(|(i, face)| {
                face.iter()
                    .all(|node| (first..first + 4).contains(node))
                    .then_some(i)
            })
            .collect();
        assert_eq!(faces.len(), 4);
        ClosedTetCavity::new(surface.clone(), &faces, cx).unwrap()
    })
}

fn bindings<'a>(
    cavities: &'a [ClosedTetCavity; 2],
    gases: &'a [ClosedGasParcel; 2],
    ambient: f64,
) -> [AdiabaticCavity<'a>; 2] {
    std::array::from_fn(|i| AdiabaticCavity {
        cavity: &cavities[i],
        reference: &gases[i],
        exterior_pressure_pa: ambient,
        controls: controls(),
    })
}

pub(super) fn problem<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    material: &'a Hyperelastic,
    dead: &'a [(usize, f64)],
) -> HyperTetProblem<'a> {
    HyperTetProblem {
        nodes_m: nodes,
        tetrahedra: tets,
        material: HyperelasticLawField::Uniform(material),
        prescribed_m: &FIXED,
        nodal_forces_n: dead,
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings {
            load_steps: 8,
            ..HyperTetSettings::default()
        },
    }
}

/// Independent entropy integral for the authored source Cv/R = 2.5 + a*T.
/// Bisection here evaluates this closed scalar formula only; production uses
/// the original NASA-9 entropy closure and complete material state admission.
fn temperature_at_ratio(slope: f64, volume_ratio: f64) -> f64 {
    let (mut lo, mut hi) = (300.0, 300.0 * volume_ratio.powf(-0.4));
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        let residual = 2.5 * (t / 300.0).ln() + slope * (t - 300.0) + volume_ratio.ln();
        if residual < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn g1_connected_cavities_keep_distinct_isentropes_in_one_loaded_equilibrium() {
    with_cx(|_, cx| {
        let (nodes, tets) = connected_shells();
        let cavities = cavities(&nodes, &tets, cx);
        let (mu, lambda, stretch) = (2.0e6, 3.0e6, 0.96_f64);
        let material =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        let piola = mu * (stretch - 1.0 / stretch) + 3.0 * lambda * stretch.ln() / stretch;
        let pressure = -piola / stretch.powi(2);
        let volume_ratio = stretch.powi(3);
        let slopes = [0.0005, 0.002];
        let temperatures = slopes.map(|slope| temperature_at_ratio(slope, volume_ratio));
        let gases = std::array::from_fn(|i| {
            gas(
                cavities[i].reference_volume_m3(),
                pressure * volume_ratio * 300.0 / temperatures[i],
                slopes[i],
                if i == 0 { 0.0 } else { 1e24 },
            )
        });
        let original = gases.clone();
        let ambient = 1.5e5;
        let bindings = bindings(&cavities, &gases, ambient);

        // The outer bipyramid volume is the sum of two positively oriented
        // tetrahedra. Their common-face contributions cancel in this oracle.
        let mut outer_gradient = vec![0.0; 3 * nodes.len()];
        for tet in [[0, 1, 2, 3], [0, 2, 1, 8]] {
            let (_, local) = tetra_volume_gradient(tet.map(|node| nodes[node]));
            for (i, &node) in tet.iter().enumerate() {
                for axis in 0..3 {
                    outer_gradient[3 * node + axis] += local[3 * i + axis];
                }
            }
        }
        let acceleration = [0.1, -0.2, -9.8];
        let field = |x: [f64; 3]| {
            Ok(GravitySample3 {
                potential_j_kg: 1e20 - dot(acceleration, x),
                acceleration_m_s2: acceleration,
                acceleration_gradient_s2: [[0.0; 3]; 3],
            })
        };
        let increment = |a: [f64; 3], b: [f64; 3]| Ok(-dot(acceleration, sub(b, a)));
        let density: Vec<_> = (0..tets.len()).map(|i| 1500.0 + 50.0 * i as f64).collect();
        let mut weight = vec![0.0; 3 * nodes.len()];
        for (tet, &rho) in tets.iter().zip(&density) {
            let x = tet.map(|node| nodes[node]);
            let quarter_mass =
                rho * dot(sub(x[1], x[0]), cross(sub(x[2], x[0]), sub(x[3], x[0]))) / 24.0;
            assert!(quarter_mass > 0.0);
            for &node in tet {
                for axis in 0..3 {
                    weight[3 * node + axis] += quarter_mass * acceleration[axis];
                }
            }
        }
        let traction = piola + ambient * stretch.powi(2);
        let dead: Vec<_> = (0..3 * nodes.len())
            .filter(|dof| !FIXED.iter().any(|&(fixed, _)| fixed == *dof))
            .map(|dof| (dof, traction * outer_gradient[dof] - weight[dof]))
            .collect();
        let problem = problem(&nodes, &tets, &material, &dead);
        let target_u: Vec<_> = nodes
            .iter()
            .flatten()
            .map(|x| (stretch - 1.0) * x)
            .collect();
        let internal = problem.internal_force(&target_u, cx).unwrap();
        let mut target_external = weight.clone();
        for &(dof, value) in &dead {
            target_external[dof] += value;
        }
        let pressure_force = cavities[0]
            .surface()
            .force_n(&target_u, TetPressure::Uniform(ambient), cx)
            .unwrap();
        for (force, pressure) in target_external.iter_mut().zip(pressure_force) {
            *force += pressure;
        }
        for source in &bindings {
            let load = source.evaluate(&target_u, cx).unwrap();
            for (force, cavity) in target_external.iter_mut().zip(load.force_n) {
                *force += cavity;
            }
        }
        for dof in 0..target_u.len() {
            if !FIXED.iter().any(|&(fixed, _)| fixed == dof) {
                near(internal[dof], target_external[dof], 2e-6);
            }
        }
        let mut forward: Option<Vec<[f64; 3]>> = None;
        for order in [[0, 1], [1, 0]] {
            let ordered = order.map(|i| bindings[i]);
            let solution = problem
                .solve_with_adiabatic_cavities(
                    &ordered,
                    Some((
                        ReferenceGravityDensity::PerElement(&density),
                        StationaryGravity3::with_potential_increment(&field, &increment),
                    )),
                    Some((cavities[0].surface(), TetPressure::Uniform(ambient))),
                    cx,
                )
                .unwrap();
            assert_eq!(solution.cavities.len(), 2);
            let solid = &solution.solid;
            for (&actual, &reference) in solid
                .displacement_m
                .iter()
                .flatten()
                .zip(nodes.iter().flatten())
            {
                near(actual, (stretch - 1.0) * reference, 5e-8);
            }
            near(solid.min_det_f, volume_ratio, 5e-8);
            if let Some(forward) = &forward {
                for (&a, &b) in solid
                    .displacement_m
                    .iter()
                    .flatten()
                    .zip(forward.iter().flatten())
                {
                    near(a, b, 5e-8);
                }
            } else {
                forward = Some(solid.displacement_m.clone());
            }
            for (slot, &source) in order.iter().enumerate() {
                let endpoint = &solution.cavities[slot];
                near(
                    endpoint.volume_m3,
                    volume_ratio * cavities[source].reference_volume_m3(),
                    4e-9,
                );
                near(endpoint.net_pressure_pa + ambient, pressure, 0.15);
                near(
                    endpoint.gas.state().state().temperature,
                    temperatures[source],
                    4e-6,
                );
                assert_eq!(
                    endpoint.gas.mass_kg().to_bits(),
                    gases[source].mass_kg().to_bits()
                );
                assert_eq!(endpoint.gas.model(), gases[source].model());
                let temperature = endpoint.gas.state().state().temperature;
                let caloric = gases[source].mass_kg()
                    * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                    / 0.028
                    * (temperature - 300.0)
                    * (2.5 + 0.5 * slopes[source] * (temperature + 300.0));
                let expected = caloric
                    + ambient * (endpoint.volume_m3 - cavities[source].reference_volume_m3());
                near(endpoint.potential_change_j, expected, 2e-8);
            }
            for &(dof, reaction) in &solid.reactions_n {
                near(reaction, traction * outer_gradient[dof] - weight[dof], 0.2);
            }
            for axis in 0..3 {
                let reaction: f64 = solid
                    .reactions_n
                    .iter()
                    .filter(|(dof, _)| dof % 3 == axis)
                    .map(|(_, force)| force)
                    .sum();
                let dead: f64 = dead
                    .iter()
                    .filter(|(dof, _)| dof % 3 == axis)
                    .map(|(_, force)| force)
                    .sum();
                let gravity: f64 = weight.iter().skip(axis).step_by(3).sum();
                near(reaction + dead + gravity, 0.0, 0.3);
            }
        }
        assert_eq!(gases, original, "all reference parcels remain immutable");
    });
}

#[test]
fn g1_heating_one_sealed_chamber_changes_the_other_through_the_connected_solid() {
    with_cx(|_, cx| {
        let (nodes, tets) = connected_shells();
        let cavities = cavities(&nodes, &tets, cx);
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2e6,
                lambda: 3e6,
            },
            4.0,
        )
        .unwrap();
        let problem = problem(&nodes, &tets, &material, &[]);
        let ambient = 1e5;
        let gases = [
            gas(cavities[0].reference_volume_m3(), ambient, 0.001, 0.0),
            gas(cavities[1].reference_volume_m3(), ambient, 0.002, 1e24),
        ];
        let baseline = problem
            .solve_with_adiabatic_cavities(&bindings(&cavities, &gases, ambient), None, None, cx)
            .unwrap();
        assert_eq!(baseline.solid.displacement_m, vec![[0.0; 3]; nodes.len()]);
        let (heated, supplied_heat) = gases[0].at_temperature(330.0).unwrap();
        assert!(supplied_heat > 0.0);
        let changed = [heated, gases[1].clone()];
        let saved = changed.clone();
        let response = problem
            .solve_with_adiabatic_cavities(&bindings(&cavities, &changed, ambient), None, None, cx)
            .unwrap();
        assert!(response.cavities[0].volume_m3 > baseline.cavities[0].volume_m3);
        assert!(
            (response.cavities[1].volume_m3 - baseline.cavities[1].volume_m3).abs() > 1e-9,
            "the second gas volume must respond through shared solid degrees of freedom"
        );
        assert!((response.cavities[1].gas.state().state().temperature - 300.0).abs() > 1e-6);
        for i in 0..2 {
            let endpoint = &response.cavities[i];
            let temperature = endpoint.gas.state().state().temperature;
            let reference_temperature = changed[i].state().state().temperature;
            let entropy_over_r = 2.5 * (temperature / reference_temperature).ln()
                + [0.001, 0.002][i] * (temperature - reference_temperature)
                + (endpoint.volume_m3 / changed[i].volume_m3()).ln();
            near(entropy_over_r, 0.0, 1e-11);
            assert_eq!(
                endpoint.gas.mass_kg().to_bits(),
                gases[i].mass_kg().to_bits()
            );
            assert_eq!(endpoint.gas.model(), gases[i].model());
        }
        assert_eq!(changed, saved);
    });
}

#[test]
fn g0_multiple_cavity_reductions_and_geometry_ownership_are_explicit() {
    with_cx(|gate, cx| {
        let (nodes, tets) = connected_shells();
        let cavities = cavities(&nodes, &tets, cx);
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2e6,
                lambda: 3e6,
            },
            4.0,
        )
        .unwrap();
        let problem = problem(&nodes, &tets, &material, &[]);
        let empty = problem
            .solve_with_adiabatic_cavities::<ClosedGasParcel>(&[], None, None, cx)
            .unwrap();
        assert_eq!(empty.solid, problem.solve(cx).unwrap());
        assert!(empty.cavities.is_empty());
        let gas = gas(cavities[0].reference_volume_m3(), 1.05e5, 0.001, 1e24);
        let source = AdiabaticCavity {
            cavity: &cavities[0],
            reference: &gas,
            exterior_pressure_pa: 1e5,
            controls: controls(),
        };
        let single = problem
            .solve_with_adiabatic_cavities(&[source], None, None, cx)
            .unwrap();
        assert_eq!(
            single.solid,
            problem
                .solve_with_adiabatic_cavity_loads(source, None, None, cx)
                .unwrap()
        );
        let u: Vec<_> = single
            .solid
            .displacement_m
            .iter()
            .flatten()
            .copied()
            .collect();
        let evaluated = source.evaluate(&u, cx).unwrap();
        assert_eq!(single.cavities.len(), 1);
        assert_eq!(single.cavities[0].gas, evaluated.gas);
        assert_eq!(single.cavities[0].force_n, evaluated.force_n);
        assert_eq!(
            single.cavities[0].potential_change_j,
            evaluated.potential_change_j
        );

        // The second owner is independently prepared, so duplicate inventory
        // must be rejected by boundary identity rather than pointer equality.
        let duplicate = ClosedTetCavity::new(
            surface(&nodes, &tets, cx),
            cavities[0].solid_face_indices(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            problem.solve_with_adiabatic_cavities(
                &[
                    source,
                    AdiabaticCavity {
                        cavity: &duplicate,
                        ..source
                    }
                ],
                None,
                None,
                cx,
            ),
            Err(HyperTetError::InvalidInput { .. })
        ));
        let translated: Vec<_> = nodes.iter().map(|x| [x[0] + 1.0, x[1], x[2]]).collect();
        let foreign = ClosedTetCavity::new(
            surface(&translated, &tets, cx),
            cavities[0].solid_face_indices(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            problem.solve_with_adiabatic_cavities(
                &[AdiabaticCavity {
                    cavity: &foreign,
                    ..source
                }],
                None,
                None,
                cx,
            ),
            Err(HyperTetError::InvalidInput { .. })
        ));
        gate.request();
        assert!(matches!(
            problem.solve_with_adiabatic_cavities(&[source], None, None, cx),
            Err(HyperTetError::Cancelled)
        ));
    });
}
