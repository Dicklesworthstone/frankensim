//! Full anisotropic source tensors must drive the original Kirchhoff/rod owner.
use super::*;
use fs_blake3::ContentHash;
use fs_couple::piezoelectric::{PiezoelectricRod, SmallSignalLimits};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_material::piezoelectric::tensor::TensorPiezoelectric;

fn material() -> TensorPiezoelectric {
    // Deliberately non-diagonal compliance and transverse coupling. The axial
    // reduction is s33=0.1, d33=0.2, epsilon33^T=2; 1/C33 and epsilon33^S
    // would give a different circuit and fail the independently derived oracle.
    let mut compliance = [[0.0; 6]; 6];
    for (i, value) in [0.125, 0.15, 0.1, 0.2, 0.25, 0.3].into_iter().enumerate() {
        compliance[i][i] = value;
    }
    compliance[0][2] = -0.02;
    compliance[2][0] = -0.02;
    compliance[1][2] = -0.03;
    compliance[2][1] = -0.03;
    let mut d = [[0.0; 6]; 3];
    d[0][5] = 0.1;
    d[1][4] = 0.1;
    d[2] = [-0.07, -0.1, 0.2, 0.0, 0.0, 0.0];
    let permittivity = [1.5, 1.7, 2.0];
    let mut claims = ClaimSet::new();
    for role in TensorPiezoelectric::strain_charge_requirements().unwrap() {
        let name = role.name();
        let value = if name == "density" {
            3.0
        } else {
            let index = name.rsplit('_').next().unwrap().as_bytes();
            let (row, col) = ((index[0] - b'1') as usize, (index[1] - b'1') as usize);
            if name.starts_with("piezoelectric_compliance") {
                compliance[row][col]
            } else if name.starts_with("piezoelectric_strain_charge") {
                d[row][col]
            } else if row == col {
                permittivity[row]
            } else {
                0.0
            }
        };
        let dims = role.dims();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored anisotropic tensor/circuit regression".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic anisotropic ceramic".into(),
            phase: "solid".into(),
            process: "full strain-charge source, no separate scalar card".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    TensorPiezoelectric::resolve_strain_charge(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        ContentHash([93; 32]),
    )
    .unwrap()
}

#[test]
fn g1_full_tensor_drives_reciprocal_motion_current_and_load_heat() {
    let tensor = material();
    for sign in [-1.0, 1.0] {
        let law = tensor.longitudinal([0.0, 0.0, sign]).unwrap();
        assert_eq!(law.material(), tensor.material());
        assert!(law.tensor_reduction().is_some());
        let identity = law.identity();
        let rod = PiezoelectricRod::new(
            law,
            2.0,
            0.5,
            SmallSignalLimits {
                max_abs_strain: 0.05,
                max_abs_electric_field_v_m: 10.0,
            },
        )
        .unwrap();
        near(rod.short_circuit_stiffness_n_m(), 2.5, 2e-14);
        near(rod.axially_clamped_capacitance_f(), 0.4, 2e-14);
        near(rod.coupling_c_m(), sign * 0.5, 2e-14);
        let (mut x, mut v, mut q) = (0.02, 0.01, sign * 0.02);
        let mut body = PiezoelectricCircuit::new(
            CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::Capacitor { farads: f64::NAN }),
                    (1, 0, Branch::Resistor { ohms: f64::NAN }),
                ],
                transformers: vec![],
            },
            vec![source(1, 2.0)],
            vec![],
            vec![q],
            0,
            rod,
            x,
            v,
            &[],
            || false,
        )
        .unwrap();
        let initial_energy = energy(x, v, q, sign * 0.2, 3.0);
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..12 {
            let expected = load_oracle(x, v, q, sign * 0.2, 3.0, 2.0, 0.03, 0.05);
            let next = body
                .step_candidate(vec![source(1, 2.0)], &[], 0.03, 0.05, solve(), || false)
                .unwrap();
            near(next.model().extension_m(), expected.0, 2e-10);
            near(next.model().velocity_m_s(), expected.1, 2e-9);
            near(next.model().charge_c(), expected.2, 2e-10);
            near(next.step().resistors[0].voltage_v, expected.3, 2e-10);
            near(next.step().heat_j, 0.05 * expected.3.powi(2) / 2.0, 2e-10);
            heat += next.step().heat_j;
            work += next.step().mechanical_supply_j;
            near(
                next.step().stored_after_j + heat - work,
                initial_energy,
                2e-10,
            );
            assert_eq!(next.model().rod().material().identity(), identity);
            assert_eq!(next.model().rod().material().material(), tensor.material());
            body = next.into_parts().0;
            (x, v, q) = (expected.0, expected.1, expected.2);
        }
        assert!(heat > 0.0);
    }
}

#[test]
fn g3_crystal_and_cut_corotation_preserves_the_native_circuit_trajectory() {
    let tensor = material();
    // A proper rotation sending original +z to world +x.
    let rotated = tensor
        .rotated(
            ContentHash([94; 32]),
            [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],
        )
        .unwrap();
    let run = |law| {
        let rod = PiezoelectricRod::new(
            law,
            2.0,
            0.5,
            SmallSignalLimits {
                max_abs_strain: 0.05,
                max_abs_electric_field_v_m: 10.0,
            },
        )
        .unwrap();
        let body = PiezoelectricCircuit::new(
            CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::Capacitor { farads: f64::NAN }),
                    (1, 0, Branch::Resistor { ohms: f64::NAN }),
                ],
                transformers: vec![],
            },
            vec![source(1, 2.0)],
            vec![],
            vec![0.02],
            0,
            rod,
            0.02,
            0.01,
            &[],
            || false,
        )
        .unwrap();
        body.step_candidate(vec![source(1, 2.0)], &[], 0.03, 0.05, solve(), || false)
            .unwrap()
    };
    let a = run(tensor.longitudinal([0.0, 0.0, 1.0]).unwrap());
    let b = run(rotated.longitudinal([1.0, 0.0, 0.0]).unwrap());
    near(a.model().extension_m(), b.model().extension_m(), 2e-12);
    near(a.model().velocity_m_s(), b.model().velocity_m_s(), 2e-12);
    near(a.model().charge_c(), b.model().charge_c(), 2e-12);
    near(a.step().heat_j, b.step().heat_j, 2e-12);
    assert_ne!(
        a.model().rod().material().identity(),
        b.model().rod().material().identity()
    );
    assert_eq!(
        a.model().rod().material().material(),
        b.model().rod().material().material()
    );
}
