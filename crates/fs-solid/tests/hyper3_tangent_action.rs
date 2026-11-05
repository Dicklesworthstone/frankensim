//! G0/G1/G3/G4: exact native material actions without global stiffness assembly.

use std::cell::Cell;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::MaterialError;
use fs_material::hyper::cards::{FiniteStrainPotential, HyperelasticMaterial};
use fs_material::hyper::{Hyperelastic, HyperelasticModel, Volumetric};
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
const TANGENT_BYTES: usize = 2 * 81 * 8;

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

fn problem<L: FiniteStrainPotential>(material: &[L; 2]) -> HyperTetProblem<'_, L> {
    HyperTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        material: HyperelasticLawField::PerElement(material),
        prescribed_m: &[(0, 0.0)],
        nodal_forces_n: &[(3, 7.0)],
        body_force_n_m3: [1.0, -2.0, 3.0],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings::default(),
    }
}

fn materials() -> [Hyperelastic; 2] {
    [
        Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 3.0,
            },
            4.0,
        )
        .unwrap(),
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 1.5,
                c01: 0.4,
                kappa: 120.0,
            },
            4.0,
        )
        .unwrap(),
    ]
}

fn displacement() -> [f64; 15] {
    // Distinct gradients on the two cells, finite shear and finite prestress.
    [
        0.0, 0.0, 0.0, 0.25, 0.0, 0.125, 0.125, -0.125, 0.0, 0.0, 0.0, 0.125, -0.125, 0.0625, 0.0,
    ]
}

fn direction() -> [f64; 15] {
    core::array::from_fn(|i| ((i * 7 + 3) % 11) as f64 * 0.0625 - 0.25)
}

fn max_abs(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0, |scale, value| scale.max(value.abs()))
}

fn close(actual: &[f64], expected: &[f64], relative: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = max_abs(expected);
    assert!(scale > 0.0);
    for (&a, &b) in actual.iter().zip(expected) {
        assert!(
            a.is_finite() && (a - b).abs() <= relative * scale,
            "{a:e} != {b:e}"
        );
    }
}

#[test]
fn heterogeneous_actions_match_native_matrix_and_independent_force_derivative() {
    with_cx(|_, cx| {
        let laws = materials();
        let p = problem(&laws);
        let prepared = p.prepare(cx).unwrap();
        let u = displacement();
        let tangent = prepared.linearized_tangent(&u, TANGENT_BYTES, cx).unwrap();
        assert_eq!(tangent.dimension(), u.len());
        let matrix = prepared.tangent_matrix(&u, cx).unwrap();
        let mut scratch = [f64::NAN; 15];
        let mut actual = [91.0; 15];
        let mut cached = [93.0; 15];
        let mut expected = [0.0; 15];
        for factor in [1e-24, 1.0, -4.0] {
            let increment = direction().map(|x| factor * x);
            prepared
                .tangent_apply(&u, &increment, &mut actual, &mut scratch, cx)
                .unwrap();
            tangent
                .apply(&increment, &mut cached, &mut scratch, cx)
                .unwrap();
            assert_eq!(actual, cached);
            matrix.spmv(&increment, &mut expected);
            close(&actual, &expected, 2e-13);
            for i in 0..3 {
                let resultant: f64 = (0..5).map(|node| actual[3 * node + i]).sum();
                assert!(resultant.abs() <= 2e-13 * max_abs(&actual));
            }
        }
        let increment = direction();
        let h = 2e-6;
        let plus = core::array::from_fn::<_, 15, _>(|i| u[i] + h * increment[i]);
        let minus = core::array::from_fn::<_, 15, _>(|i| u[i] - h * increment[i]);
        let force_plus = prepared.internal_force(&plus, cx).unwrap();
        let force_minus = prepared.internal_force(&minus, cx).unwrap();
        let derivative: Vec<_> = force_plus
            .iter()
            .zip(force_minus)
            .map(|(a, b)| (a - b) / (2.0 * h))
            .collect();
        tangent
            .apply(&increment, &mut actual, &mut scratch, cx)
            .unwrap();
        close(&actual, &derivative, 2e-8);
        // The material-only operator retains constrained rows and is independent
        // of prescribed dead forces, which have no displacement derivative.
        assert!(actual[0].abs() > 1e-3);
    });
}

#[test]
fn gent_chain_stiffening_drives_original_solid_force_and_tangent() {
    with_cx(|_, cx| {
        let parameters = [(2.4e5, 6.0), (4.1e5, 9.0)];
        let laws = parameters.map(|(mu, locking_invariant)| {
            Hyperelastic::new(
                HyperelasticModel::Gent {
                    mu,
                    locking_invariant,
                    volumetric: Volumetric::Quadratic { kappa: 4.0e6 },
                },
                8.0,
            )
            .unwrap()
        });
        let prepared = problem(&laws).prepare(cx).unwrap();
        // Both original cells undergo F=I+gamma e_x tensor e_y, hence J=1.
        // Node 2 is the only nonzero shear-coordinate shape value. Its x force
        // and tangent action equal the derivatives of sum_e V_e W_e(gamma).
        let shear_direction =
            core::array::from_fn::<_, 15, _>(|d| if d % 3 == 0 { NODES[d / 3][1] } else { 0.0 });
        let mut action = [0.0; 15];
        let mut scratch = [0.0; 15];
        for gamma in [0.2, 1.2, 2.35] {
            let u = shear_direction.map(|value| gamma * value);
            let (energy, force, stiffness) =
                parameters
                    .iter()
                    .fold((0.0, 0.0, 0.0), |(energy, force, stiffness), &(mu, jm)| {
                        let ratio = gamma * gamma / jm;
                        (
                            energy - mu * jm * (-ratio).ln_1p() / 12.0,
                            force + mu * gamma / (6.0 * (1.0 - ratio)),
                            stiffness + mu * (1.0 + ratio) / (6.0 * (1.0 - ratio).powi(2)),
                        )
                    });
            let native_force = prepared.internal_force(&u, cx).unwrap();
            assert!((native_force[6] - force).abs() <= 4e-11 * force.abs());
            assert!(
                (prepared.stored_energy_j(&u, cx).unwrap() - energy).abs() <= 4e-11 * energy.abs()
            );
            let tangent = prepared.linearized_tangent(&u, TANGENT_BYTES, cx).unwrap();
            tangent
                .apply(&shear_direction, &mut action, &mut scratch, cx)
                .unwrap();
            assert!((action[6] - stiffness).abs() <= 4e-11 * stiffness.abs());
            let mut assembled = [0.0; 15];
            prepared
                .tangent_matrix(&u, cx)
                .unwrap()
                .spmv(&shear_direction, &mut assembled);
            close(&action, &assembled, 2e-13);
        }
        let locked = shear_direction.map(|value| 2.5 * value);
        let before = action;
        assert!(matches!(
            prepared.tangent_apply(&locked, &shear_direction, &mut action, &mut scratch, cx),
            Err(HyperTetError::MaterialRefused { element: 0, .. })
        ));
        assert_eq!(action, before);
        assert!(prepared.internal_force(&locked, cx).is_err());
    });
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    [-v[1], v[0], v[2]]
}

#[test]
fn current_frame_covariance_and_exact_translation_null_action() {
    with_cx(|_, cx| {
        let laws = materials();
        let prepared = problem(&laws).prepare(cx).unwrap();
        let u = displacement();
        let du = direction();
        let mut output = [0.0; 15];
        let mut scratch = [0.0; 15];
        prepared
            .tangent_apply(&u, &du, &mut output, &mut scratch, cx)
            .unwrap();
        let rotated_u = core::array::from_fn::<_, 15, _>(|d| {
            let node = d / 3;
            let current = core::array::from_fn(|i| NODES[node][i] + u[3 * node + i]);
            rotate(current)[d % 3] - NODES[node][d % 3] + [4.0, -2.0, 1.0][d % 3]
        });
        let rotated_du = core::array::from_fn::<_, 15, _>(|d| {
            rotate(core::array::from_fn(|i| du[3 * (d / 3) + i]))[d % 3]
        });
        let expected = core::array::from_fn::<_, 15, _>(|d| {
            rotate(core::array::from_fn(|i| output[3 * (d / 3) + i]))[d % 3]
        });
        let operator = prepared
            .linearized_tangent(&rotated_u, TANGENT_BYTES, cx)
            .unwrap();
        operator
            .apply(&rotated_du, &mut output, &mut scratch, cx)
            .unwrap();
        close(&output, &expected, 2e-13);
        let translation = core::array::from_fn::<_, 15, _>(|d| {
            [1_073_741_824.0, -536_870_912.0, 268_435_456.0][d % 3]
        });
        operator
            .apply(&translation, &mut output, &mut scratch, cx)
            .unwrap();
        assert_eq!(output, [0.0; 15]);
        operator
            .apply_transpose(&translation, &mut output, &mut scratch, cx)
            .unwrap();
        assert_eq!(output, [0.0; 15]);
    });
}

/// Independent homogeneous thermal-reference pullback W(F)=s^3 W_el(F/s).
/// Its elastic-domain bound is intentionally narrower than the TOTAL stretch.
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
fn original_reference_adapter_retains_thermal_tangent_and_elastic_domain() {
    with_cx(|_, cx| {
        let laws = [1.25, 1.25].map(|stretch| ReferenceDilation {
            law: Hyperelastic::new(
                HyperelasticModel::NeoHookean {
                    mu: 2.0,
                    lambda: 3.0,
                },
                1.05,
            )
            .unwrap(),
            stretch,
        });
        let prepared = problem(&laws).prepare(cx).unwrap();
        let u = core::array::from_fn::<_, 15, _>(|d| 0.25 * NODES[d / 3][d % 3]);
        let f = [1.25, 0.0, 0.0, 0.0, 1.25, 0.0, 0.0, 0.0, 1.25];
        assert!(laws[0].law.tangent(&f).is_err());
        let mut increment = [0.0; 15];
        increment[3] = 0.5;
        let mut output = [0.0; 15];
        let mut scratch = [0.0; 15];
        prepared
            .tangent_apply(&u, &increment, &mut output, &mut scratch, cx)
            .unwrap();
        assert!((output[3] - 1.25 * (3.0 + 2.0 * 2.0) * 0.5 / 6.0).abs() < 2e-13);
        let operator = prepared.linearized_tangent(&u, TANGENT_BYTES, cx).unwrap();
        let mut cached = [0.0; 15];
        operator
            .apply(&increment, &mut cached, &mut scratch, cx)
            .unwrap();
        assert_eq!(output, cached);
        let matrix = prepared.tangent_matrix(&u, cx).unwrap();
        matrix.spmv(&increment, &mut cached);
        close(&output, &cached, 2e-13);
    });
}

/// A deliberately nonsymmetric numerical callback probes the actual transpose;
/// it is an adapter test, not an admissible physical elastic-energy model.
struct TangentProbe<'a> {
    calls: Cell<usize>,
    cancel: Option<&'a CancelGate>,
    bad: bool,
}
impl FiniteStrainPotential for TangentProbe<'_> {
    fn energy_density(&self, _: &[f64; 9]) -> f64 {
        panic!("tangent must not query energy")
    }
    fn first_piola(&self, _: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        panic!("tangent must not query stress")
    }
    fn consistent_tangent(&self, _: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        self.calls.set(self.calls.get() + 1);
        if let Some(gate) = self.cancel {
            gate.request();
        }
        let mut matrix =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 2.0 } else { 0.0 }));
        matrix[0][4] = if self.bad { f64::NAN } else { 0.75 };
        Ok(matrix)
    }
}

#[test]
fn frozen_actions_reuse_original_point_hessians_and_apply_actual_transpose() {
    with_cx(|_, cx| {
        let laws = core::array::from_fn(|_| TangentProbe {
            calls: Cell::new(0),
            cancel: None,
            bad: false,
        });
        let prepared = problem(&laws).prepare(cx).unwrap();
        let u = displacement();
        assert!(
            prepared
                .linearized_tangent(&u, TANGENT_BYTES - 1, cx)
                .is_err()
        );
        assert_eq!(laws[0].calls.get() + laws[1].calls.get(), 0);
        let operator = prepared.linearized_tangent(&u, TANGENT_BYTES, cx).unwrap();
        assert_eq!(laws[0].calls.get() + laws[1].calls.get(), 2);
        let mut output = [0.0; 15];
        let mut transposed = [0.0; 15];
        let mut scratch = [0.0; 15];
        let x = direction();
        let y = core::array::from_fn::<_, 15, _>(|i| x[14 - i]);
        for _ in 0..4 {
            operator.apply(&x, &mut output, &mut scratch, cx).unwrap();
            operator
                .apply_transpose(&y, &mut transposed, &mut scratch, cx)
                .unwrap();
        }
        assert_eq!(laws[0].calls.get() + laws[1].calls.get(), 2);
        let xy: f64 = output.iter().zip(y).map(|(a, b)| a * b).sum();
        let yx: f64 = transposed.iter().zip(x).map(|(a, b)| a * b).sum();
        assert!((xy - yx).abs() <= 2e-13 * xy.abs().max(1.0));
        let matrix = prepared.tangent_matrix(&u, cx).unwrap();
        let expected: Vec<_> = (0..15)
            .map(|i| (0..15).map(|j| matrix.get(j, i) * y[j]).sum())
            .collect();
        close(&transposed, &expected, 2e-13);
        operator.apply(&y, &mut output, &mut scratch, cx).unwrap();
        assert!(
            output
                .iter()
                .zip(transposed)
                .any(|(a, b)| (a - b).abs() > 1e-3)
        );
    });
}

#[test]
fn input_domain_overflow_and_cancellation_never_publish_partial_actions() {
    with_cx(|gate, cx| {
        let laws = materials();
        let prepared = problem(&laws).prepare(cx).unwrap();
        let u = displacement();
        let operator = prepared.linearized_tangent(&u, TANGENT_BYTES, cx).unwrap();
        let mut output = [17.0; 15];
        let mut scratch = [0.0; 15];
        assert!(
            prepared
                .tangent_apply(&u[..14], &direction(), &mut output, &mut scratch, cx)
                .is_err()
        );
        assert!(
            operator
                .apply(&direction()[..14], &mut output, &mut scratch, cx)
                .is_err()
        );
        assert!(
            operator
                .apply(&[f64::NAN; 15], &mut output, &mut scratch, cx)
                .is_err()
        );
        assert!(
            operator
                .apply(&direction(), &mut output[..14], &mut scratch, cx)
                .is_err()
        );
        assert!(
            operator
                .apply(&direction(), &mut output, &mut scratch[..14], cx)
                .is_err()
        );
        assert_eq!(output, [17.0; 15]);
        let mut inverted_second = u;
        inverted_second[12] = 2.0;
        assert!(matches!(
            prepared.tangent_apply(
                &inverted_second,
                &direction(),
                &mut output,
                &mut scratch,
                cx
            ),
            Err(HyperTetError::MaterialRefused { element: 1, .. })
        ));
        assert!(
            prepared
                .linearized_tangent(&inverted_second, TANGENT_BYTES, cx)
                .is_err()
        );
        assert_eq!(output, [17.0; 15]);
        let mut overflow = [f64::MAX; 15];
        overflow[0] = -f64::MAX;
        assert!(
            operator
                .apply(&overflow, &mut output, &mut scratch, cx)
                .is_err()
        );
        assert_eq!(output, [17.0; 15]);
        gate.request();
        assert!(matches!(
            operator.apply(&direction(), &mut output, &mut scratch, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            operator.apply_transpose(&direction(), &mut output, &mut scratch, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            prepared.tangent_apply(&u, &direction(), &mut output, &mut scratch, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            prepared.linearized_tangent(&u, TANGENT_BYTES, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert_eq!(output, [17.0; 15]);
    });
    for bad in [false, true] {
        with_cx(|gate, cx| {
            let laws = [
                TangentProbe {
                    calls: Cell::new(0),
                    cancel: (!bad).then_some(gate),
                    bad: false,
                },
                TangentProbe {
                    calls: Cell::new(0),
                    cancel: None,
                    bad,
                },
            ];
            let prepared = problem(&laws).prepare(cx).unwrap();
            let mut output = [19.0; 15];
            let result = prepared.tangent_apply(
                &displacement(),
                &direction(),
                &mut output,
                &mut [0.0; 15],
                cx,
            );
            if bad {
                assert!(matches!(
                    result,
                    Err(HyperTetError::MaterialRefused { element: 1, .. })
                ));
                assert!(
                    prepared
                        .linearized_tangent(&displacement(), TANGENT_BYTES, cx)
                        .is_err()
                );
            } else {
                assert!(matches!(result, Err(HyperTetError::Cancelled)));
                assert_eq!(laws[1].calls.get(), 0);
            }
            assert_eq!(output, [19.0; 15]);
        });
    }
}
