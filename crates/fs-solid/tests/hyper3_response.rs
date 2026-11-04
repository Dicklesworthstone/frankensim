//! G0/G1/G3/G4: original finite-strain fields, work conjugacy and refusals.

use std::cell::Cell;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::MaterialError;
use fs_material::hyper::cards::{FiniteStrainPotential, HyperelasticMaterial};
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, Volumetric};
use fs_solid::{
    HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField, TetAssemblyBudget,
};

const NODES: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [-1.0, 0.0, 0.0],
];
const TETS: [[usize; 4]; 2] = [[0, 1, 2, 3], [0, 4, 3, 2]];
const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

fn problem<'a, L: FiniteStrainPotential>(
    nodes: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    material: HyperelasticLawField<'a, L>,
) -> HyperTetProblem<'a, L> {
    HyperTetProblem {
        nodes_m: nodes,
        tetrahedra,
        material,
        prescribed_m: &[],
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings::default(),
    }
}

fn nh(mu: f64, lambda: f64, limit: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, limit).unwrap()
}

fn affine(nodes: &[[f64; 3]], f: [f64; 9], translation: [f64; 3]) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|x| {
            core::array::from_fn::<_, 3, _>(|i| {
                translation[i]
                    + (0..3)
                        .map(|j| (f[3 * i + j] - if i == j { 1.0 } else { 0.0 }) * x[j])
                        .sum::<f64>()
            })
        })
        .collect()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= 5e-11 * scale.max(1.0),
        "{actual:e} != {expected:e}, scale {scale:e}"
    );
}

fn product(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

fn transpose(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| a[3 * (ij % 3) + ij / 3])
}

#[test]
fn heterogeneous_fields_match_closed_form_stress_and_original_nodal_forces() {
    with_cx(|_, cx| {
        let parameters = [(2.0, 3.0), (7.0, 5.0)];
        let materials = parameters.map(|(mu, lambda)| nh(mu, lambda, 4.0));
        let p = problem(&NODES, &TETS, HyperelasticLawField::PerElement(&materials));
        let stretch = [1.2_f64, 0.8, 1.1];
        let f = [
            stretch[0], 0.0, 0.0, 0.0, stretch[1], 0.0, 0.0, 0.0, stretch[2],
        ];
        let u = affine(&NODES, f, [7.0, -11.0, 5.0]);
        let prepared = p.prepare(cx).unwrap();
        let fields = prepared.element_responses(&u, cx).unwrap();
        assert_eq!(fields, p.element_responses(&u, cx).unwrap());
        assert_eq!(fields.len(), 2);
        let j = stretch.iter().product::<f64>();
        let log_j = j.ln();
        let gradients = [
            [
                [-1.0, -1.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            [
                [1.0, -1.0, -1.0],
                [-1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0],
            ],
        ];
        let mut reconstructed_force = vec![0.0; 3 * NODES.len()];
        for (element, response) in fields.iter().enumerate() {
            let (mu, lambda) = parameters[element];
            assert_eq!(response.element, element);
            assert_eq!(response.nodes, TETS[element]);
            near(response.jacobian, j, 1.0);
            near(response.reference_volume_m3, 1.0 / 6.0, 1.0);
            near(response.current_volume_m3, j / 6.0, 1.0);
            for i in 0..3 {
                for k in 0..3 {
                    let piola = if i == k {
                        mu * stretch[i] + (lambda * log_j - mu) / stretch[i]
                    } else {
                        0.0
                    };
                    let cauchy = if i == k {
                        (mu * stretch[i] * stretch[i] + lambda * log_j - mu) / j
                    } else {
                        0.0
                    };
                    near(response.first_piola_pa[3 * i + k], piola, mu + lambda);
                    near(response.cauchy_pa[3 * i + k], cauchy, mu + lambda);
                    near(
                        response.green_lagrange_strain[3 * i + k],
                        if i == k {
                            0.5 * (stretch[i] * stretch[i] - 1.0)
                        } else {
                            0.0
                        },
                        1.0,
                    );
                }
            }
            let potential = 0.5 * mu * (stretch.iter().map(|s| s * s).sum::<f64>() - 3.0)
                - mu * log_j
                + 0.5 * lambda * log_j * log_j;
            near(response.potential_density_j_m3, potential, mu + lambda);
            near(response.potential_j, potential / 6.0, mu + lambda);
            for (a, &node) in response.nodes.iter().enumerate() {
                for i in 0..3 {
                    reconstructed_force[3 * node + i] += response.reference_volume_m3
                        * (0..3)
                            .map(|k| response.first_piola_pa[3 * i + k] * gradients[element][a][k])
                            .sum::<f64>();
                }
            }
        }
        // Shared interface nodes collect BOTH original material contributions.
        // Recovery must not smooth their distinct cell stresses into one field.
        assert!((fields[0].cauchy_pa[0] - fields[1].cauchy_pa[0]).abs() > 1.0);
        let native_force = prepared.internal_force(&u, cx).unwrap();
        for (actual, expected) in reconstructed_force.iter().zip(native_force) {
            near(*actual, expected, 20.0);
        }
        near(
            fields.iter().map(|response| response.potential_j).sum(),
            prepared.stored_energy_j(&u, cx).unwrap(),
            20.0,
        );
    });
}

#[test]
fn anisotropic_spatial_stress_is_objective_and_preserves_current_area_traction() {
    with_cx(|_, cx| {
        let material = Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 5.0,
                k1: 12.0,
                k2: 1.1,
                dispersion: 0.05,
                fibres: vec![[0.8, 0.0, 0.6], [0.0, 0.6, 0.8]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 80.0 },
            },
            4.0,
        )
        .unwrap();
        let p = problem(
            &NODES[..4],
            &TETS[..1],
            HyperelasticLawField::Uniform(&material),
        );
        let prepared = p.prepare(cx).unwrap();
        let f = [1.2, 0.1, 0.0, 0.0, 0.9, 0.08, 0.04, 0.0, 1.1];
        let rotation = [0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let u = affine(&NODES[..4], f, [0.0; 3]);
        let rotated_u = affine(&NODES[..4], product(&rotation, &f), [6.0, -2.0, 1.0]);
        let before = prepared.element_responses(&u, cx).unwrap().remove(0);
        let after = prepared
            .element_responses(&rotated_u, cx)
            .unwrap()
            .remove(0);
        let rotated_piola = product(&rotation, &before.first_piola_pa);
        let rotated_cauchy = product(
            &product(&rotation, &before.cauchy_pa),
            &transpose(&rotation),
        );
        for ij in 0..9 {
            near(after.first_piola_pa[ij], rotated_piola[ij], 100.0);
            near(after.cauchy_pa[ij], rotated_cauchy[ij], 100.0);
            near(
                after.green_lagrange_strain[ij],
                before.green_lagrange_strain[ij],
                1.0,
            );
            near(
                before.cauchy_pa[ij],
                before.cauchy_pa[3 * (ij % 3) + ij / 3],
                100.0,
            );
        }
        near(after.potential_j, before.potential_j, 100.0);
        near(after.current_volume_m3, before.current_volume_m3, 1.0);

        // Independent Nanson map: sigma*(cof(F)*A0) = P*A0, with arbitrary
        // oriented reference area A0. Piola and Cauchy must use different areas.
        let cof = [
            f[4] * f[8] - f[5] * f[7],
            f[5] * f[6] - f[3] * f[8],
            f[3] * f[7] - f[4] * f[6],
            f[2] * f[7] - f[1] * f[8],
            f[0] * f[8] - f[2] * f[6],
            f[1] * f[6] - f[0] * f[7],
            f[1] * f[5] - f[2] * f[4],
            f[2] * f[3] - f[0] * f[5],
            f[0] * f[4] - f[1] * f[3],
        ];
        let reference_area = [0.7, -0.2, 0.4];
        let current_area: [f64; 3] =
            core::array::from_fn(|i| (0..3).map(|j| cof[3 * i + j] * reference_area[j]).sum());
        for i in 0..3 {
            let current_force: f64 = (0..3)
                .map(|j| before.cauchy_pa[3 * i + j] * current_area[j])
                .sum();
            let reference_force: f64 = (0..3)
                .map(|j| before.first_piola_pa[3 * i + j] * reference_area[j])
                .sum();
            near(current_force, reference_force, 100.0);
        }

        let rigid_u = affine(&NODES[..4], rotation, [6.0, -2.0, 1.0]);
        let rigid = prepared.element_responses(&rigid_u, cx).unwrap().remove(0);
        for value in rigid.green_lagrange_strain {
            near(value, 0.0, 1.0);
        }
        for value in rigid.cauchy_pa {
            near(value, 0.0, 100.0);
        }
        near(rigid.potential_j, 0.0, 100.0);
    });
}

struct ReferenceDilation {
    law: Hyperelastic,
    stretch: f64,
}

impl HyperelasticMaterial for ReferenceDilation {
    fn hyperelastic(&self) -> &Hyperelastic {
        &self.law
    }

    fn energy_density(&self, f: &[f64; 9]) -> f64 {
        self.stretch.powi(3) * self.law.energy(&f.map(|value| value / self.stretch))
    }

    fn first_piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        Ok(self
            .law
            .piola(&f.map(|value| value / self.stretch))?
            .map(|value| value * self.stretch.powi(2)))
    }

    fn consistent_tangent(&self, f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        Ok(self
            .law
            .tangent(&f.map(|value| value / self.stretch))?
            .map(|row| row.map(|value| value * self.stretch)))
    }
}

#[test]
fn recovery_retains_reference_deformation_adapters_and_reports_total_strain() {
    with_cx(|_, cx| {
        let material = ReferenceDilation {
            law: nh(2.0, 3.0, 1.05),
            stretch: 1.3,
        };
        let f = IDENTITY.map(|value| value * material.stretch);
        assert!(material.law.piola(&f).is_err());
        let p = problem(
            &NODES[..4],
            &TETS[..1],
            HyperelasticLawField::Uniform(&material),
        );
        let u = affine(&NODES[..4], f, [0.0; 3]);
        let response = p.element_responses(&u, cx).unwrap().remove(0);
        for ij in 0..9 {
            near(response.deformation_gradient[ij], f[ij], 1.0);
            near(response.first_piola_pa[ij], 0.0, 1.0);
            near(response.cauchy_pa[ij], 0.0, 1.0);
            near(
                response.green_lagrange_strain[ij],
                if ij / 3 == ij % 3 {
                    0.5 * (material.stretch.powi(2) - 1.0)
                } else {
                    0.0
                },
                1.0,
            );
        }
        near(response.potential_j, 0.0, 1.0);
        near(
            response.current_volume_m3,
            material.stretch.powi(3) / 6.0,
            1.0,
        );
    });
}

/// A frozen potential with an arbitrary energy datum and a combined response.
/// The separate methods deliberately refuse: recovery must not duplicate a
/// constitutive update, request its tangent or relabel its potential as heat.
struct CombinedPotential {
    law: Hyperelastic,
    calls: Cell<usize>,
    datum: f64,
    bad_stress: bool,
}

impl CombinedPotential {
    fn new(datum: f64, bad_stress: bool) -> Self {
        Self {
            law: nh(2.0, 3.0, 4.0),
            calls: Cell::new(0),
            datum,
            bad_stress,
        }
    }
}

impl FiniteStrainPotential for CombinedPotential {
    fn energy_density(&self, _f: &[f64; 9]) -> f64 {
        panic!("recovery must use the original combined response")
    }

    fn first_piola(&self, _f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        panic!("recovery must use the original combined response")
    }

    fn consistent_tangent(&self, _f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        panic!("stress recovery does not need a constitutive tangent")
    }

    fn stress_and_potential(&self, f: &[f64; 9]) -> Result<([f64; 9], f64), MaterialError> {
        self.calls.set(self.calls.get() + 1);
        let mut piola = self.law.piola(f)?;
        if self.bad_stress {
            piola[0] = f64::NAN;
        }
        Ok((piola, self.law.energy(f) + self.datum))
    }
}

#[test]
fn sparse_recovery_observes_only_original_nodes_and_keeps_incremental_potential() {
    with_cx(|_, cx| {
        let materials = [
            CombinedPotential::new(-3.0, false),
            CombinedPotential::new(-7.0, false),
        ];
        let p = problem(&NODES, &TETS, HyperelasticLawField::PerElement(&materials));
        let prepared = p.prepare(cx).unwrap();
        let mut sampled = Vec::new();
        let response = prepared
            .element_response_with(
                1,
                |node| {
                    sampled.push(node);
                    [0.0; 3]
                },
                cx,
            )
            .unwrap();
        assert_eq!(sampled, TETS[1]);
        assert_eq!(response.element, 1);
        assert_eq!(response.nodes, TETS[1]);
        assert_eq!(materials[0].calls.get(), 0);
        assert_eq!(materials[1].calls.get(), 1);
        near(response.potential_density_j_m3, -7.0, 7.0);
        near(response.potential_j, -7.0 / 6.0, 7.0);
        let complete = prepared.element_responses(&[0.0; 15], cx).unwrap();
        assert_eq!(response, complete[1]);
        assert_eq!(materials[0].calls.get(), 1);
        assert_eq!(materials[1].calls.get(), 2);
        near(complete[0].potential_density_j_m3, -3.0, 3.0);
    });
}

#[test]
fn recovery_refuses_invalid_states_nonfinite_material_fields_and_cancellation() {
    with_cx(|gate, cx| {
        let material = nh(2.0, 3.0, 4.0);
        let p = problem(
            &NODES[..4],
            &TETS[..1],
            HyperelasticLawField::Uniform(&material),
        );
        let prepared = p.prepare(cx).unwrap();
        assert!(prepared.element_responses(&[0.0; 11], cx).is_err());
        assert!(prepared.element_responses(&[f64::NAN; 12], cx).is_err());
        for stretch in [-1.0, 5.0] {
            let mut f = IDENTITY;
            f[0] = stretch;
            let u = affine(&NODES[..4], f, [0.0; 3]);
            assert!(matches!(
                prepared.element_responses(&u, cx),
                Err(HyperTetError::MaterialRefused { element: 0, .. })
            ));
        }
        assert!(
            prepared
                .element_response_with(1, |_| panic!("missing element must not sample"), cx,)
                .is_err()
        );
        assert!(
            prepared
                .element_response_with(0, |_| [f64::NAN, 0.0, 0.0], cx)
                .is_err()
        );
        for (datum, bad_stress) in [(f64::NAN, false), (0.0, true)] {
            let materials = [
                CombinedPotential::new(0.0, false),
                CombinedPotential::new(datum, bad_stress),
            ];
            let bad = problem(&NODES, &TETS, HyperelasticLawField::PerElement(&materials));
            assert!(matches!(
                bad.element_responses(&[0.0; 15], cx),
                Err(HyperTetError::MaterialRefused { element: 1, .. })
            ));
        }
        gate.request();
        assert!(matches!(
            prepared.element_responses(&[0.0; 12], cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            prepared.element_response_with(
                0,
                |_| panic!("cancelled observation must not sample"),
                cx,
            ),
            Err(HyperTetError::Cancelled)
        ));
    });
}
