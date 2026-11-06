//! G0/G3: material-axis covariance of native Hill plasticity and its energy.

use fs_material::plastic::energy::{EnergyTolerances, IncrementEnergy};
use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::{Hill48, Hill48Coefficients, HillResponse};
use fs_material::{IsotropicElastic, SmallStrainLaw, Voigt};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
// Proper rotation from the unit quaternion (1,2,3,4)/sqrt(30).
const Q: [[f64; 3]; 3] = [
    [-2.0 / 3.0, 2.0 / 15.0, 11.0 / 15.0],
    [2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0],
    [1.0 / 3.0, 14.0 / 15.0, 2.0 / 15.0],
];
const FIRST: Voigt = [0.004, -0.001, -0.0008, 0.0015, 0.0003, -0.0004];
const SECOND: Voigt = [0.006, -0.0015, 0.0007, -0.0011, 0.0013, 0.0002];
const TOLERANCE: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-9,
    yield_pa: 1e-8,
};

fn law() -> Hill48 {
    Hill48::new(
        IsotropicElastic::new(1000.0, 0.3, 0.05).unwrap(),
        1.0,
        60.0,
        Hill48Coefficients {
            f: 0.35,
            g: 0.8,
            h: 0.45,
            l: 1.5,
            m: 1.1,
            n: 2.0,
        },
    )
    .unwrap()
}

// Independent full-tensor rotation; do not use the production Voigt transform.
fn rotated(v: Voigt, q: [[f64; 3]; 3]) -> Voigt {
    let tensor = [[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]];
    let value = |i: usize, j: usize| -> f64 {
        (0..3)
            .flat_map(|a| (0..3).map(move |b| q[i][a] * tensor[a][b] * q[j][b]))
            .sum()
    };
    [
        value(0, 0),
        value(1, 1),
        value(2, 2),
        value(0, 1),
        value(1, 2),
        value(2, 0),
    ]
}

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= allowance,
        "actual {actual:.16e}, expected {expected:.16e}, allowance {allowance:e}"
    );
}

fn energies(e: &IncrementEnergy) -> [f64; 10] {
    [
        e.initial.elastic_j_m3,
        e.initial.hardening_j_m3,
        e.initial.kinematic_hardening_j_m3,
        e.final_state.elastic_j_m3,
        e.final_state.hardening_j_m3,
        e.final_state.kinematic_hardening_j_m3,
        e.mechanical_work_j_m3,
        e.physical_dissipation_j_m3,
        e.numerical_dissipation_j_m3,
        e.residual_j_m3,
    ]
}

fn response_bits(r: &HillResponse) -> Vec<u64> {
    r.stress_pa
        .iter()
        .chain(r.tangent_pa.iter().flatten())
        .chain(&r.state.plastic_strain)
        .chain([&r.state.alpha, &r.plastic_increment, &r.yield_residual_pa])
        .map(|v| v.to_bits())
        .collect()
}

#[test]
fn rotated_load_history_preserves_native_plastic_history_and_work_partition() {
    let native = law();
    let oriented = OrientedHill48::new(native.clone(), Q).unwrap();
    assert_eq!(oriented.material_to_spatial(), &Q);
    let (mut local_state, mut spatial_state) = (native.initial_state(), oriented.initial_state());
    let (mut old_local, mut old_spatial) = ([0.0; 6], [0.0; 6]);
    for local_strain in [
        FIRST,
        SECOND,
        [0.002, -0.0006, -0.001, 0.0008, -0.0004, 0.0014],
    ] {
        let spatial_strain = rotated(local_strain, Q);
        let a = native
            .advance_with_energy(&old_local, &local_state, &local_strain, TOLERANCE)
            .unwrap();
        let b = oriented
            .advance_with_energy(&old_spatial, &spatial_state, &spatial_strain, TOLERANCE)
            .unwrap();
        assert_eq!(b.strain, spatial_strain);
        let expected_stress = rotated(a.response.stress_pa, Q);
        let expected_plastic = rotated(a.response.state.plastic_strain, Q);
        let observed_plastic = oriented.plastic_strain_spatial(&b.response.state);
        for j in 0..6 {
            near(b.response.stress_pa[j], expected_stress[j], 2e-9);
            near(
                b.response.state.plastic_strain[j],
                a.response.state.plastic_strain[j],
                2e-12,
            );
            near(observed_plastic[j], expected_plastic[j], 2e-12);
        }
        near(b.response.state.alpha, a.response.state.alpha, 2e-12);
        near(
            b.response.plastic_increment,
            a.response.plastic_increment,
            2e-12,
        );
        // Native solvers re-evaluate accepted endpoints before the next step;
        // rotation roundoff must not cause an unrepresentable plastic return.
        let repeated = oriented
            .response(&spatial_strain, &b.response.state)
            .unwrap();
        near(repeated.plastic_increment, 0.0, 2e-12);
        near(repeated.state.alpha, b.response.state.alpha, 2e-12);
        for (actual, expected) in energies(&b.energy).into_iter().zip(energies(&a.energy)) {
            near(actual, expected, 2e-11);
        }
        let spatial_work: f64 = (0..6)
            .map(|j| {
                let weight = if j < 3 { 1.0 } else { 2.0 };
                weight * b.response.stress_pa[j] * (spatial_strain[j] - old_spatial[j])
            })
            .sum();
        near(b.energy.mechanical_work_j_m3, spatial_work, 2e-11);
        near(
            oriented
                .stored_energy(&spatial_strain, &b.response.state)
                .unwrap()
                .total_j_m3(),
            native
                .stored_energy(&local_strain, &a.response.state)
                .unwrap()
                .total_j_m3(),
            2e-11,
        );
        assert!(b.energy.physical_dissipation_j_m3 >= 0.0);
        local_state = a.response.state;
        spatial_state = b.response.state;
        old_local = local_strain;
        old_spatial = spatial_strain;
    }
    assert!(local_state.alpha > 0.0);
}

#[test]
fn off_axis_plastic_tangent_is_work_symmetric_and_matches_all_six_derivatives() {
    let native = law();
    let old = native
        .response(&FIRST, &native.initial_state())
        .unwrap()
        .state;
    let oriented = OrientedHill48::new(native, Q).unwrap();
    let strain = rotated(SECOND, Q);
    let response = oriented.response(&strain, &old).unwrap();
    assert!(response.plastic_increment > 0.0);
    let mut worst = 0.0_f64;
    let mut scale = 0.0_f64;
    let h = 1e-8;
    for j in 0..6 {
        let (mut plus, mut minus) = (strain, strain);
        plus[j] += h;
        minus[j] -= h;
        let plus = oriented.response(&plus, &old).unwrap();
        let minus = oriented.response(&minus, &old).unwrap();
        for i in 0..6 {
            let finite_difference = (plus.stress_pa[i] - minus.stress_pa[i]) / (2.0 * h);
            worst = worst.max((response.tangent_pa[i][j] - finite_difference).abs());
            scale = scale.max(response.tangent_pa[i][j].abs());
            let wi = if i < 3 { 1.0 } else { 2.0 };
            let wj = if j < 3 { 1.0 } else { 2.0 };
            near(
                wi * response.tangent_pa[i][j],
                wj * response.tangent_pa[j][i],
                2e-8,
            );
        }
    }
    assert!(
        worst / scale < 2e-7,
        "off-axis tangent relative error {}",
        worst / scale
    );
}

#[test]
fn identity_orientation_preserves_every_response_and_energy_bit() {
    let native = law();
    let oriented = OrientedHill48::new(native.clone(), IDENTITY).unwrap();
    let mut state = native.initial_state();
    let mut old = [0.0; 6];
    for strain in [FIRST, SECOND, [0.0; 6]] {
        let a = native
            .advance_with_energy(&old, &state, &strain, TOLERANCE)
            .unwrap();
        let b = oriented
            .advance_with_energy(&old, &state, &strain, TOLERANCE)
            .unwrap();
        assert_eq!(response_bits(&a.response), response_bits(&b.response));
        assert_eq!(
            energies(&a.energy).map(f64::to_bits),
            energies(&b.energy).map(f64::to_bits)
        );
        assert_eq!(
            oriented
                .plastic_strain_spatial(&b.response.state)
                .map(f64::to_bits),
            a.response.state.plastic_strain.map(f64::to_bits)
        );
        state = a.response.state;
        old = strain;
    }
}

#[test]
fn invalid_material_axes_are_refused_without_normalizing_them() {
    for (i, j, value) in [
        (0, 0, -1.0),
        (1, 1, 1.01),
        (0, 1, 0.01),
        (2, 2, 0.0),
        (2, 0, f64::NAN),
        (1, 2, f64::INFINITY),
    ] {
        let mut axes = IDENTITY;
        axes[i][j] = value;
        assert!(OrientedHill48::new(law(), axes).is_err());
    }
    assert!(
        OrientedHill48::new(law(), [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]).is_ok()
    );
}
