use super::*;
use crate::elastic::acoustic::interface::surface::{
    ElasticSurfaceCondition, ElasticSurfaceImpedance, ElasticSurfaceImpedanceTensor,
    IsotropicElasticSurface,
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const POLARIZATIONS: [[f64; 3]; 3] = [
    [4.0 / 5.0, 3.0 / 5.0, 0.0],
    [-36.0 / 65.0, 48.0 / 65.0, 5.0 / 13.0],
    [3.0 / 13.0, -4.0 / 13.0, 12.0 / 13.0],
];

// Analytically constructed SPD Mandel tensor. Its normal Christoffel matrix
// has eigenpairs rho*c^2 and POLARIZATIONS, independent of the native eigensolver.
fn anisotropic(ceiling: f64) -> ElasticAcousticState {
    let diagonal = [
        20e9, 18e9, 40_211_502_958.579_88, 12e9,
        29_985_276_213.017_75, 21_591_717_869.822_487,
    ];
    let mut c = [[0.0; 6]; 6];
    for i in 0..6 { c[i][i] = diagonal[i]; }
    for (i, j, v) in [
        (2, 4, -11_307_817_337.162_804),
        (2, 5, 8_480_863_002.872_103),
        (4, 5, -14_388_957_159.763_313),
    ] {
        c[i][j] = v;
        c[j][i] = v;
    }
    ElasticAcousticState::try_new_mandel(c, 2400.0, ContentHash([71; 32]), ceiling).unwrap()
}

fn isotropic() -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        2500.0, 2500.0 * (4000.0_f64.powi(2) - (4.0 / 3.0) * 2000.0_f64.powi(2)),
        2500.0 * 2000.0_f64.powi(2), ContentHash([73; 32]), 0.02,
    ).unwrap()
}

fn support() -> ElasticSurfaceImpedanceTensor {
    let mut z = [[C64::ZERO; 3]; 3];
    z[0][0] = C64::new(2_664_359.861_591_695_5, 3_557_093.425_605_536_3);
    z[1][1] = C64::new(3e6, -1e6);
    z[2][2] = C64::new(4_335_640.138_408_304, 2_442_906.574_394_463_7);
    z[0][2] = C64::new(1_245_674.740_484_429, -830_449.826_989_619_4);
    z[2][0] = z[0][2];
    ElasticSurfaceImpedanceTensor::new(z).unwrap()
}

fn incoming_cartesian() -> [C64; 3] {
    [C64::from_re(0.001), C64::new(-0.0002, 0.0001), C64::new(0.0003, -0.0004)]
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:?} != {expected:?}");
}

fn scalar_close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:?} != {expected:?}");
}

#[test]
fn original_christoffel_characteristic_keeps_tensor_modes_strain_and_ideal_phase() {
    let state = anisotropic(0.02);
    let owner = NormalElasticCharacteristic::new(&state, IDENTITY).unwrap();
    assert!(core::ptr::eq(owner.state(), &state));
    assert_eq!(owner.surface_to_material(), IDENTITY);
    assert!(!state.is_explicitly_isotropic());
    let speeds = [1500.0, 2600.0, 4300.0];
    for m in 0..3 {
        let wave = owner.incoming_mode(m).unwrap();
        let reverse = owner.reflected_mode(m).unwrap();
        scalar_close(owner.directional_waves().modes[m].phase_speed_m_s, speeds[m], 1e-7);
        scalar_close(owner.mode_impedances_pa_s_m()[m], 2400.0 * speeds[m], 2e-4);
        for i in 0..3 {
            close(wave.velocity[i], C64::from_re(POLARIZATIONS[m][i]), 2e-10);
            assert_eq!(wave.velocity[i], reverse.velocity[i]);
            assert_eq!(wave.traction_pa_s_m[i], -reverse.traction_pa_s_m[i]);
            close(wave.traction_pa_s_m[i], wave.velocity[i].scale(-2400.0 * speeds[m]), 2e-3);
        }
        for i in 0..6 {
            assert_eq!(wave.strain_mandel_s_m[i], -reverse.strain_mandel_s_m[i]);
        }
        let stress: [f64; 6] = core::array::from_fn(|i| {
            (0..6).map(|j| state.stiffness_mandel_pa()[i][j] * wave.strain_mandel_s_m[j].re).sum()
        });
        let traction = [
            stress[5] * core::f64::consts::FRAC_1_SQRT_2,
            stress[4] * core::f64::consts::FRAC_1_SQRT_2,
            stress[2],
        ];
        for i in 0..3 { close(wave.traction_pa_s_m[i], C64::from_re(traction[i]), 2e-7); }
        scalar_close(owner.mode_energy_velocities_surface()[m][2], speeds[m], 2e-7);
    }
    assert!(owner.incoming_mode(3).is_none());
    let d = owner.impedance_pa_s_m();
    let y = owner.admittance_m_pa_s();
    for i in 0..3 { for j in 0..3 {
        let product: f64 = (0..3).map(|k| d[i][k] * y[k][j]).sum();
        scalar_close(product, if i == j { 1.0 } else { 0.0 }, 3e-10);
    } }
    let a = [C64::new(0.001, -0.0002), C64::new(-0.0003, 0.0007), C64::from_re(0.0005)];
    for condition in [ElasticSurfaceCondition::TractionFree, ElasticSurfaceCondition::Clamped] {
        let prepared = owner.prepare_surface(condition).unwrap();
        let response = prepared.drive(a).unwrap();
        for m in 0..3 {
            let sign = if condition == ElasticSurfaceCondition::Clamped { -1.0 } else { 1.0 };
            assert_eq!(response.reflected_velocity_m_s[m], a[m].scale(sign));
            assert_eq!(prepared.velocity_coefficients(m).unwrap()[m], C64::from_re(sign));
        }
        assert_eq!(response.absorbed_normal_power_w_m2, 0.0);
        scalar_close(response.total_incoming_power_w_m2, response.total_reflected_power_w_m2, 1e-12);
        let constrained = if condition == ElasticSurfaceCondition::Clamped {
            response.boundary_velocity_m_s
        } else { response.boundary_traction_pa };
        assert_eq!(constrained, [C64::ZERO; 3]);
    }
}

// Independent 70-digit physical Cartesian solve of (D+Z)v=2D a. D and Z do
// not commute: diagonalizing the support or reversing matrix order is wrong.
#[test]
fn anisotropic_noncommuting_support_matches_independent_physical_oracle() {
    let state = anisotropic(0.02);
    let owner = NormalElasticCharacteristic::new(&state, IDENTITY).unwrap();
    let prepared = owner.prepare_surface(ElasticSurfaceCondition::TensorImpedance(support())).unwrap();
    let response = prepared.drive_cartesian(incoming_cartesian()).unwrap();
    let expected_in = [
        C64::new(0.00068, 0.00006),
        C64::new(-0.0005861538461538462, -0.00008),
        C64::new(0.0005692307692307692, -0.0004),
    ];
    let expected_reflected = [
        C64::new(-0.0000988057541760178, -0.0004162447230261288),
        C64::new(-0.00023694744358727365, 0.0001787377395451169),
        C64::new(0.00011364304399667953, -0.00016850700780363323),
    ];
    let expected_wall = [
        C64::new(0.0010784131447990634, -0.00047087522058519085),
        C64::new(-0.00046922711669211416, 0.0000340924224033731),
        C64::new(0.00031376763923259893, -0.00048679964583984725),
    ];
    for i in 0..3 {
        close(response.incoming_velocity_m_s[i], expected_in[i], 2e-13);
        close(response.reflected_velocity_m_s[i], expected_reflected[i], 2e-13);
        close(response.boundary_velocity_m_s[i], expected_wall[i], 2e-13);
    }
    scalar_close(response.total_incoming_power_w_m2, 4.428288284023669, 3e-9);
    scalar_close(response.total_reflected_power_w_m2, 0.8174410699969257, 3e-9);
    scalar_close(response.absorbed_normal_power_w_m2, 3.610847214026743, 3e-9);
    scalar_close(response.boundary_normal_flux_w_m2, response.absorbed_normal_power_w_m2, 3e-9);
    let load = support();
    let expected_traction = load.traction(response.boundary_velocity_m_s);
    for i in 0..3 { close(response.boundary_traction_pa[i], expected_traction[i], 2e-8); }
    assert_eq!(response.material_state_identity, state.material_state_identity());
}

#[test]
fn common_frame_rotation_and_degenerate_isotropic_basis_preserve_cartesian_response() {
    let q = [[0.6, 0.0, 0.8], [0.0, 1.0, 0.0], [-0.8, 0.0, 0.6]];
    let state = anisotropic(0.02);
    let rotated = ElasticAcousticState::try_new_mandel(
        crate::tensor::rotate_mandel(*state.stiffness_mandel_pa(), q),
        state.density_kg_m3(), ContentHash([74; 32]), 0.02,
    ).unwrap();
    let condition = ElasticSurfaceCondition::TensorImpedance(support());
    let a = NormalElasticCharacteristic::new(&state, IDENTITY).unwrap()
        .prepare_surface(condition).unwrap().drive_cartesian(incoming_cartesian()).unwrap();
    let b = NormalElasticCharacteristic::new(&rotated, q).unwrap()
        .prepare_surface(condition).unwrap().drive_cartesian(incoming_cartesian()).unwrap();
    for i in 0..3 {
        close(a.boundary_velocity_m_s[i], b.boundary_velocity_m_s[i], 3e-13);
        close(a.boundary_traction_pa[i], b.boundary_traction_pa[i], 3e-6);
    }
    scalar_close(a.absorbed_normal_power_w_m2, b.absorbed_normal_power_w_m2, 3e-9);
    let transform = crate::tensor::mandel_rotation(q);
    for i in 0..6 {
        let strain = (0..6).fold(C64::ZERO, |sum, j| sum + a.boundary_strain_mandel[j].scale(transform[i][j]));
        close(b.boundary_strain_mandel[i], strain, 2e-15);
    }
    // The old P polarization reverses on reflection. Comparing physical
    // fields instead of falsely equating modal indices resolves that sign.
    let isotropic = isotropic();
    let input = incoming_cartesian();
    let old = IsotropicElasticSurface::new(&isotropic, condition).unwrap()
        .prepare_incidence(ElasticMode::P, 1.0).unwrap()
        .drive([input[2], input[0], input[1]]).unwrap();
    for frame in [IDENTITY, q] {
        let response = NormalElasticCharacteristic::new(&isotropic, frame).unwrap()
            .prepare_surface(condition).unwrap().drive_cartesian(input).unwrap();
        for i in 0..3 {
            close(response.boundary_velocity_m_s[i], old.boundary_velocity_m_s[i], 3e-13);
            close(response.boundary_traction_pa[i], old.boundary_traction_pa[i], 3e-6);
        }
        scalar_close(response.absorbed_normal_power_w_m2, old.absorbed_normal_power_w_m2, 3e-9);
    }
}

#[test]
fn separate_reflection_and_wall_maps_keep_matching_and_stiff_support_work() {
    let state = isotropic();
    let owner = NormalElasticCharacteristic::new(&state, IDENTITY).unwrap();
    let z = owner.mode_impedances_pa_s_m()[2];
    let a = C64::from_re(0.001);
    for resistance in [z * (1.0 + 1e-9), 1e28] {
        let support = ElasticSurfaceImpedance::new(
            C64::from_re(resistance), C64::from_re(resistance),
        ).unwrap();
        let response = owner.prepare_surface(ElasticSurfaceCondition::Impedance(support))
            .unwrap().drive([C64::ZERO, C64::ZERO, a]).unwrap();
        let expected_reflection = a.scale((z - resistance) / (z + resistance));
        let expected_velocity = a.scale(2.0 * (z / resistance) / (1.0 + z / resistance));
        close(response.reflected_velocity_m_s[2], expected_reflection, expected_reflection.abs() * 3e-14);
        close(response.boundary_velocity_m_s[2], expected_velocity, expected_velocity.abs() * 3e-14);
        let expected_loss = 0.5 * resistance * expected_velocity.norm_sq();
        assert!(response.absorbed_normal_power_w_m2 > 0.0);
        scalar_close(response.absorbed_normal_power_w_m2, expected_loss, expected_loss * 5e-14);
    }
    let reactive = ElasticSurfaceImpedance::new(C64::new(0.0, 1e8), C64::new(0.0, -2e8)).unwrap();
    let response = owner.prepare_surface(ElasticSurfaceCondition::Impedance(reactive)).unwrap()
        .drive([C64::from_re(0.0002), C64::ZERO, a]).unwrap();
    assert_eq!(response.absorbed_normal_power_w_m2, 0.0);
    scalar_close(response.total_reflected_power_w_m2, response.total_incoming_power_w_m2, 1e-10);
}

#[test]
fn normal_owner_enforces_original_halfspace_strain_frame_and_cancellation() {
    let state = anisotropic(1e-8);
    let owner = NormalElasticCharacteristic::new(&state, IDENTITY).unwrap();
    // Free-surface strain cancels at z=0, but the standing wave away from the
    // boundary still exceeds this original source's domain.
    let free = owner.prepare_surface(ElasticSurfaceCondition::TractionFree).unwrap();
    assert!(free.drive([C64::from_re(0.001), C64::ZERO, C64::ZERO]).is_err());
    assert!(free.drive([C64::ZERO; 3]).is_ok());
    assert!(free.drive([C64::from_re(f64::NAN), C64::ZERO, C64::ZERO]).is_err());
    for value in [1.01, -1.0, f64::NAN] {
        let mut frame = IDENTITY;
        frame[0][0] = value;
        assert!(NormalElasticCharacteristic::new(&state, frame).is_err());
    }
    assert!(NormalElasticCharacteristic::new_cancellable(&state, IDENTITY, &mut || true).is_err());
    assert!(owner.prepare_surface_cancellable(
        ElasticSurfaceCondition::TensorImpedance(support()), &mut || true,
    ).is_err());
}

fn second_material(ceiling: f64) -> ElasticAcousticState {
    let p = [[15.0 / 17.0, 0.0, -8.0 / 17.0], [0.0, 1.0, 0.0], [8.0 / 17.0, 0.0, 15.0 / 17.0]];
    let speeds = [1800.0_f64, 2300.0, 5000.0];
    let gamma: [[f64; 3]; 3] = core::array::from_fn(|i| core::array::from_fn(|j| {
        (0..3).map(|m| 3000.0 * speeds[m].powi(2) * p[m][i] * p[m][j]).sum()
    }));
    let mut c = [[0.0; 6]; 6];
    c[0][0] = 20e9;
    c[1][1] = 18e9;
    c[3][3] = 12e9;
    let axes = [2, 1, 0];
    let indices = [2, 4, 5];
    let scales = [1.0, core::f64::consts::SQRT_2, core::f64::consts::SQRT_2];
    for i in 0..3 { for j in i..3 {
        let value = scales[i] * scales[j] * gamma[axes[i]][axes[j]];
        c[indices[i]][indices[j]] = value;
        c[indices[j]][indices[i]] = value;
    } }
    ElasticAcousticState::try_new_mandel(c, 3000.0, ContentHash([72; 32]), ceiling).unwrap()
}

// Independent Cartesian oracle: v=(D1+D2)^-1*2(D1*a+D2*c), b1=v-a,b2=v-c.
#[test]
fn welded_two_material_normal_modes_match_two_sided_oracle_and_keep_both_domains() {
    let first = anisotropic(0.02);
    let second = second_material(0.02);
    let a = NormalElasticCharacteristic::new(&first, IDENTITY).unwrap();
    let b = NormalElasticCharacteristic::new(&second, IDENTITY).unwrap();
    let prepared = PreparedNormalElasticInterface::new(&a, &b).unwrap();
    let input = [incoming_cartesian(), [
        C64::new(-0.0001, 0.0002), C64::from_re(0.0003), C64::new(-0.0002, -0.0001),
    ]];
    let response = prepared.drive_cartesian(input).unwrap();
    let expected = [
        C64::new(0.0006245028712234452, 0.0002973839570102374),
        C64::new(-0.00008115353355547229, 0.00015901067122488372),
        C64::new(-0.00004609916141606349, -0.00045350415516028606),
    ];
    for s in 0..2 {
        let owner = prepared.characteristic(if s == 0 { ElasticSide::First } else { ElasticSide::Second });
        let outgoing = owner.modal_velocity_to_surface(response.outgoing_velocity_m_s[s]);
        for i in 0..3 {
            close(response.boundary_velocity_m_s[s][i], expected[i], 3e-13);
            close(outgoing[i], expected[i] - input[s][i], 3e-13);
        }
    }
    scalar_close(response.total_incoming_power_w_m2, 5.248788284023669, 4e-9);
    scalar_close(response.total_outgoing_power_w_m2, 5.248788284023669, 4e-9);
    for i in 0..3 {
        close(response.boundary_traction_pa[0][i], response.boundary_traction_pa[1][i], 3e-6);
    }
    assert_eq!(response.material_state_identities, [first.material_state_identity(), second.material_state_identity()]);
    let transparent = PreparedNormalElasticInterface::new(&a, &a).unwrap();
    let modal = [a.surface_velocity_to_modal(input[0]), a.surface_velocity_to_modal(input[1])];
    let unchanged = transparent.drive(modal).unwrap();
    assert_eq!(unchanged.outgoing_velocity_m_s, [modal[1], modal[0]]);
    let limited = second_material(1e-9);
    let limited = NormalElasticCharacteristic::new(&limited, IDENTITY).unwrap();
    assert!(PreparedNormalElasticInterface::new(&a, &limited).unwrap().drive_cartesian(input).is_err());
    assert!(PreparedNormalElasticInterface::new_cancellable(&a, &b, &mut || true).is_err());
}
