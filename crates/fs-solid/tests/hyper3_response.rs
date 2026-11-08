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

fn affine_velocity(nodes: &[[f64; 3]], rate: &[f64; 9], translation: [f64; 3]) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|x| {
            core::array::from_fn::<_, 3, _>(|i| {
                (0..3).fold(translation[i], |v, j| rate[3 * i + j].mul_add(x[j], v))
            })
        })
        .collect()
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

#[test]
fn native_deformation_retains_summed_sub_ulp_diagonal_strain() {
    with_cx(|_, cx| {
        // All three nonanchor shape gradients have x component one. Each
        // x-displacement contribution rounds away if added separately to I,
        // although their sum produces a representable stretch and stress.
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, -1.0, -1.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let material = nh(1.0, 0.0, 2.0);
        let p = problem(&nodes, &TETS[..1], HyperelasticLawField::Uniform(&material));
        let term = 0.5 * f64::EPSILON;
        let u = [
            0.0, 0.0, 0.0, term, 0.0, 0.0, term, 0.0, 0.0, term, 0.0, 0.0,
        ];
        let prepared = p.prepare(cx).unwrap();
        let response = prepared.element_responses(&u, cx).unwrap().remove(0);
        let stretch = 1.0 + 3.0 * term;
        assert!(stretch > 1.0);
        assert_eq!(
            response.deformation_gradient,
            [stretch, term, term, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
        );
        let extension = stretch - 1.0;
        let expected_strain = extension * (1.0 + 0.5 * extension);
        assert_eq!(response.green_lagrange_strain[0], expected_strain);

        // Independent Neo-Hookean P00 = mu*(F00 - F^-T_00), lambda = 0.
        // Use a relative gate: an ordinary unit-scale tolerance would pass
        // the original zero-strain/zero-axial-stress defect.
        let expected_stress = stretch - stretch.recip();
        assert!((response.first_piola_pa[0] / expected_stress - 1.0).abs() < 16.0 * f64::EPSILON);
        let internal = prepared.internal_force(&u, cx).unwrap();
        assert!((internal[3] / (expected_stress / 6.0) - 1.0).abs() < 16.0 * f64::EPSILON);
        assert_eq!(response, p.element_responses(&u, cx).unwrap().remove(0));
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

#[test]
fn native_motion_rates_match_affine_kinematics_and_assembled_internal_work() {
    with_cx(|_, cx| {
        let materials = [nh(2.0, 3.0, 4.0), nh(7.0, 5.0, 4.0)];
        let p = problem(&NODES, &TETS, HyperelasticLawField::PerElement(&materials));
        let prepared = p.prepare(cx).unwrap();
        let stretches = [1.25, 0.8, 1.1];
        let rates = [0.125, -0.16, 0.033];
        let f = core::array::from_fn(|ij| {
            if ij / 3 == ij % 3 {
                stretches[ij / 3]
            } else {
                0.0
            }
        });
        let rate = core::array::from_fn(|ij| if ij / 3 == ij % 3 { rates[ij / 3] } else { 0.0 });
        let u = affine(&NODES, f, [7.0, -11.0, 5.0]);
        let velocity = affine_velocity(&NODES, &rate, [3.0, -2.0, 1.0]);
        let responses = prepared
            .element_motion_responses(&u, &velocity, cx)
            .unwrap();
        assert_eq!(
            responses,
            p.element_motion_responses(&u, &velocity, cx).unwrap()
        );
        let j: f64 = stretches.iter().product();
        let dilation: f64 = (0..3).map(|i| rates[i] / stretches[i]).sum();
        for response in &responses {
            for ij in 0..9 {
                let (i, k) = (ij / 3, ij % 3);
                let spatial = if i == k { rates[i] / stretches[i] } else { 0.0 };
                let material = if i == k { rates[i] * stretches[i] } else { 0.0 };
                near(response.deformation_gradient_rate_per_s[ij], rate[ij], 1.0);
                near(response.spatial_velocity_gradient_per_s[ij], spatial, 1.0);
                near(response.spatial_strain_rate_per_s[ij], spatial, 1.0);
                near(response.green_lagrange_strain_rate_per_s[ij], material, 1.0);
            }
            near(response.jacobian_rate_per_s, j * dilation, 1.0);
            near(response.current_volume_rate_m3_s, j * dilation / 6.0, 1.0);
            near(
                response.reference_stress_power_w,
                response.spatial_stress_power_w,
                10.0,
            );
        }
        // This independently exercises the original FE force assembly, rather
        // than only comparing two contractions of the recovered stress.
        let force = prepared.internal_force(&u, cx).unwrap();
        let native_power: f64 = force.iter().zip(&velocity).map(|(f, v)| f * v).sum();
        near(
            responses.iter().map(|r| r.reference_stress_power_w).sum(),
            native_power,
            10.0,
        );
        near(
            responses.iter().map(|r| r.spatial_stress_power_w).sum(),
            native_power,
            10.0,
        );
    });
}

#[test]
fn motion_rates_are_objective_under_time_dependent_superposed_rotation() {
    with_cx(|_, cx| {
        let material = nh(2.0, 3.0, 4.0);
        let p = problem(
            &NODES[..4],
            &TETS[..1],
            HyperelasticLawField::Uniform(&material),
        );
        let prepared = p.prepare(cx).unwrap();
        let f = [1.25, 0.125, 0.0, 0.0625, 0.875, 0.125, 0.0, 0.03125, 1.125];
        let rate = [
            0.25, -0.375, 0.125, 0.0625, 0.125, -0.0625, -0.125, 0.0625, 0.1875,
        ];
        let q = [0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let spin = [0.0, -0.5, 0.25, 0.5, 0.0, -0.125, -0.25, 0.125, 0.0];
        let rotated_f = product(&q, &f);
        let rotated_rate = product(&q, &rate);
        let frame_rate = product(&spin, &rotated_f);
        let total_rate = core::array::from_fn(|ij| rotated_rate[ij] + frame_rate[ij]);
        let observe = |f, rate| {
            prepared
                .element_motion_responses(
                    &affine(&NODES[..4], f, [0.0; 3]),
                    &affine_velocity(&NODES[..4], &rate, [0.0; 3]),
                    cx,
                )
                .unwrap()
                .remove(0)
        };
        let original = observe(f, rate);
        let rotated = observe(rotated_f, total_rate);
        let expected_d = product(
            &product(&q, &original.spatial_strain_rate_per_s),
            &transpose(&q),
        );
        let rotated_l = product(
            &product(&q, &original.spatial_velocity_gradient_per_s),
            &transpose(&q),
        );
        for ij in 0..9 {
            near(rotated.spatial_strain_rate_per_s[ij], expected_d[ij], 1.0);
            near(
                rotated.spatial_velocity_gradient_per_s[ij],
                rotated_l[ij] + spin[ij],
                1.0,
            );
            near(
                rotated.green_lagrange_strain_rate_per_s[ij],
                original.green_lagrange_strain_rate_per_s[ij],
                1.0,
            );
        }
        near(
            rotated.current_volume_rate_m3_s,
            original.current_volume_rate_m3_s,
            1.0,
        );
        near(
            rotated.reference_stress_power_w,
            original.reference_stress_power_w,
            10.0,
        );
        near(
            rotated.spatial_stress_power_w,
            original.spatial_stress_power_w,
            10.0,
        );
        let rigid = observe(rotated_f, frame_rate);
        for value in rigid
            .spatial_strain_rate_per_s
            .iter()
            .chain(&rigid.green_lagrange_strain_rate_per_s)
        {
            near(*value, 0.0, 1.0);
        }
        near(rigid.jacobian_rate_per_s, 0.0, 1.0);
        near(rigid.reference_stress_power_w, 0.0, 10.0);
        near(rigid.spatial_stress_power_w, 0.0, 10.0);
    });
}

#[test]
fn motion_recovery_preserves_tiny_rates_and_common_motion_invariance() {
    with_cx(|_, cx| {
        let material = nh(2.0, 3.0, 4.0);
        let p = problem(
            &NODES[..4],
            &TETS[..1],
            HyperelasticLawField::Uniform(&material),
        );
        let prepared = p.prepare(cx).unwrap();
        let tiny = 2.0_f64.powi(-90);
        let mut rate = [0.0; 9];
        rate[0] = tiny;
        let small = prepared
            .element_motion_responses(
                &[0.0; 12],
                &affine_velocity(&NODES[..4], &rate, [0.0; 3]),
                cx,
            )
            .unwrap()
            .remove(0);
        assert_eq!(small.deformation_gradient_rate_per_s[0], tiny);
        assert_eq!(small.spatial_strain_rate_per_s[0], tiny);
        assert_eq!(small.green_lagrange_strain_rate_per_s[0], tiny);
        assert_eq!(small.jacobian_rate_per_s, tiny);

        let f = [1.25, 0.125, 0.0, 0.0, 1.125, 0.0625, 0.0, 0.0, 0.875];
        let rate = [0.125, -0.25, 0.0, 0.0625, 0.125, 0.0, 0.0, 0.0, -0.0625];
        let observe = |translation, velocity_translation| {
            prepared
                .element_motion_responses(
                    &affine(&NODES[..4], f, translation),
                    &affine_velocity(&NODES[..4], &rate, velocity_translation),
                    cx,
                )
                .unwrap()
                .remove(0)
        };
        // Dyadic physical increments remain exactly represented at these
        // translations; the observer must not manufacture strain or power.
        let position = 2.0_f64.powi(40);
        let velocity = 2.0_f64.powi(35);
        assert_eq!(
            observe([0.0; 3], [0.0; 3]),
            observe(
                [position, -position, position],
                [-velocity, velocity, velocity]
            ),
        );
    });
}

#[test]
fn sparse_motion_observation_is_single_pass_and_refuses_before_material_evaluation() {
    with_cx(|gate, cx| {
        let materials = [
            CombinedPotential::new(-3.0, false),
            CombinedPotential::new(-7.0, false),
        ];
        let p = problem(&NODES, &TETS, HyperelasticLawField::PerElement(&materials));
        let prepared = p.prepare(cx).unwrap();
        let u = affine(
            &NODES,
            [1.125, 0.125, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            [0.0; 3],
        );
        let velocity = affine_velocity(&NODES, &[0.125; 9], [1.0, 2.0, 3.0]);
        let mut sampled = Vec::new();
        let sparse = prepared
            .element_motion_response_with(
                1,
                |node| {
                    sampled.push(node);
                    (
                        core::array::from_fn(|i| u[3 * node + i]),
                        core::array::from_fn(|i| velocity[3 * node + i]),
                    )
                },
                cx,
            )
            .unwrap();
        assert_eq!(sampled, TETS[1]);
        assert_eq!(materials[0].calls.get(), 0);
        assert_eq!(materials[1].calls.get(), 1);
        assert_eq!(
            sparse,
            prepared
                .element_motion_responses(&u, &velocity, cx)
                .unwrap()[1]
        );
        assert_eq!(materials[0].calls.get(), 1);
        assert_eq!(materials[1].calls.get(), 2);
        assert!(
            prepared
                .element_motion_responses(&u, &[0.0; 14], cx)
                .is_err()
        );
        assert!(
            prepared
                .element_motion_responses(&u, &[f64::NAN; 15], cx)
                .is_err()
        );
        assert!(
            prepared
                .element_motion_response_with(2, |_| panic!("missing element"), cx)
                .is_err()
        );
        assert!(
            prepared
                .element_motion_response_with(0, |_| ([0.0; 3], [f64::INFINITY; 3]), cx)
                .is_err()
        );
        let inverted = affine(
            &NODES,
            [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            [0.0; 3],
        );
        assert!(
            prepared
                .element_motion_responses(&inverted, &velocity, cx)
                .is_err()
        );
        assert_eq!(materials[0].calls.get(), 1);
        assert_eq!(materials[1].calls.get(), 2);
        assert!(matches!(
            prepared.element_motion_response_with(
                0,
                |_| {
                    gate.request();
                    ([0.0; 3], [0.0; 3])
                },
                cx
            ),
            Err(HyperTetError::Cancelled)
        ));
        assert_eq!(materials[0].calls.get(), 1);
        assert!(matches!(
            prepared.element_motion_response_with(0, |_| panic!("cancelled"), cx),
            Err(HyperTetError::Cancelled)
        ));
    });
}

#[test]
fn motion_power_retains_the_full_original_stress_and_spin() {
    // Deliberately nonobjective potential: recovery must expose its full work,
    // not make it look objective by discarding the nonsymmetric stress part.
    struct Nonobjective;
    impl FiniteStrainPotential for Nonobjective {
        fn energy_density(&self, f: &[f64; 9]) -> f64 {
            f[1] - f[3]
        }
        fn first_piola(&self, _: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
            Ok([0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        }
        fn consistent_tangent(&self, _: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
            panic!("motion recovery must not request a tangent")
        }
    }
    with_cx(|_, cx| {
        let material = Nonobjective;
        let p = problem(
            &NODES[..4],
            &TETS[..1],
            HyperelasticLawField::Uniform(&material),
        );
        let spin = [0.0, 0.5, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0];
        let response = p
            .element_motion_responses(
                &[0.0; 12],
                &affine_velocity(&NODES[..4], &spin, [0.0; 3]),
                cx,
            )
            .unwrap()
            .remove(0);
        assert_eq!(response.spatial_strain_rate_per_s, [0.0; 9]);
        assert_eq!(response.green_lagrange_strain_rate_per_s, [0.0; 9]);
        near(response.reference_stress_power_w, 1.0 / 6.0, 1.0);
        near(response.spatial_stress_power_w, 1.0 / 6.0, 1.0);
    });
}
