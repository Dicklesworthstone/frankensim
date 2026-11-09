//! Original anisotropic modes drive the same support storage and causal engine.

#[path = "normal/coating.rs"]
mod coating;

use super::{close, solid};
use fs_blake3::ContentHash;
use fs_couple::elastic_reflection::{
    NormalOscillatorSupportedSurface, OrientedOscillatorSupport,
    OrientedOscillatorSupportedSurface, PreparedTransientElasticReflection,
};
use fs_ga::{So3, So3Tangent, Vec3};
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::interface::ElasticMode;
use fs_material::elastic::acoustic::interface::normal::NormalElasticCharacteristic;
use fs_math::c64::C64;
use fs_phs::WallPin;

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn anisotropic() -> ElasticAcousticState {
    // Synthetic SPD Mandel law. Its normal Christoffel matrix has speeds
    // 1500,2600,4300 and rational mixed normal/shear polarizations, rather than
    // a rotated isotropic approximation. The independent reference uses this
    // original tensor and a full Cartesian support solve at 70-digit precision.
    let stiffness = [
        [20e9, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 18e9, 0.0, 0.0, 0.0, 0.0],
        [
            0.0,
            0.0,
            40_211_502_958.57988,
            0.0,
            -11_307_817_337.162804,
            8_480_863_002.872103,
        ],
        [0.0, 0.0, 0.0, 12e9, 0.0, 0.0],
        [
            0.0,
            0.0,
            -11_307_817_337.162804,
            0.0,
            29_985_276_213.01775,
            -14_388_957_159.763313,
        ],
        [
            0.0,
            0.0,
            8_480_863_002.872103,
            0.0,
            -14_388_957_159.763313,
            21_591_717_869.822487,
        ],
    ];
    ElasticAcousticState::try_new_mandel(stiffness, 2400.0, ContentHash([49; 32]), 0.002).unwrap()
}

fn support() -> OrientedOscillatorSupport {
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.0, 2.0 * 0.25_f64.atan(), 0.0))).unwrap();
    let pins = core::array::from_fn(|i| WallPin {
        surface_density: 1000.0,
        stiffness_per_area: [12e9, 2e9, 8e9][i],
        resistance: [2e6, 3e6, 5e6][i],
    });
    OrientedOscillatorSupport::new(pins, rotation)
}

fn cc(a: C64, b: C64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:?} != {b:?}");
}

fn modal(characteristic: &NormalElasticCharacteristic<'_>, cartesian: [C64; 3]) -> [C64; 3] {
    characteristic.mode_polarizations_surface().map(|p| {
        (0..3).fold(C64::ZERO, |value, axis| {
            value + cartesian[axis].scale(p[axis])
        })
    })
}

#[test]
fn g1_anisotropic_supported_reflection_matches_independent_cartesian_work() {
    let material = anisotropic();
    let characteristic = NormalElasticCharacteristic::new(&material, IDENTITY).unwrap();
    let support = support();
    let owner =
        NormalOscillatorSupportedSurface::prepare(&characteristic, support, 2000.0).unwrap();
    let incoming = [
        C64::from_re(0.001),
        C64::new(-0.0002, 0.0001),
        C64::new(0.0003, -0.0004),
    ];
    let response = owner.drive(modal(&characteristic, incoming)).unwrap();
    assert!(core::ptr::eq(
        owner.reflection().characteristic().state(),
        &material
    ));
    let velocity = [
        C64::new(0.0010784131447990634, -0.00047087522058519085),
        C64::new(-0.0004692271166921142, 0.00003409242240337310),
        C64::new(0.00031376763923259893, -0.00048679964583984725),
    ];
    let traction = [
        C64::new(-4534.817589546807, -1714.4729656909046),
        C64::new(1373.5889276729694, -571.5043839022335),
        C64::new(-3501.9033948394313, 2826.208432877362),
    ];
    for axis in 0..3 {
        cc(
            response.reflection.boundary_velocity_m_s[axis],
            velocity[axis],
            5e-15,
        );
        cc(
            response.reflection.boundary_traction_pa[axis],
            traction[axis],
            5e-8,
        );
    }
    close(
        response.reflection.total_incoming_power_w_m2,
        4.428288284023669,
        3e-11,
    );
    close(
        response.reflection.total_reflected_power_w_m2,
        0.8174410699969257,
        3e-11,
    );
    close(
        response.reflection.absorbed_normal_power_w_m2,
        3.610847214026743,
        3e-11,
    );
    let principal = [
        (velocity[0].scale(15.0) - velocity[2].scale(8.0)).scale(1.0 / 17.0),
        velocity[1],
        (velocity[0].scale(8.0) + velocity[2].scale(15.0)).scale(1.0 / 17.0),
    ];
    for axis in 0..3 {
        let pin = support.principal_supports()[axis];
        cc(
            response.principal_velocity_m_s[axis],
            principal[axis],
            5e-15,
        );
        cc(
            response.principal_displacement_m[axis] * C64::new(0.0, -2000.0),
            principal[axis],
            5e-15,
        );
        close(
            response.mean_kinetic_energy_j_m2[axis],
            0.25 * pin.surface_density * principal[axis].norm_sq(),
            3e-14,
        );
        close(
            response.mean_elastic_energy_j_m2[axis],
            0.25 * pin.stiffness_per_area * principal[axis].norm_sq() / 2000.0_f64.powi(2),
            3e-14,
        );
    }
}

#[test]
fn g1_anisotropic_causal_modes_match_original_physical_generator() {
    let material = anisotropic();
    let characteristic = NormalElasticCharacteristic::new(&material, IDENTITY).unwrap();
    let owner = PreparedTransientElasticReflection::prepare_normal_characteristic(
        &characteristic,
        support(),
        0.0002,
    )
    .unwrap();
    let mut state = owner
        .state(0.0, [1e-8, -2e-8, 0.5e-8], [3e-5, -2e-5, 1e-5])
        .unwrap();
    let mut workspace = owner.workspace();
    let input = modal(
        &characteristic,
        [
            C64::from_re(0.001),
            C64::from_re(-0.0002),
            C64::from_re(0.0003),
        ],
    )
    .map(|v| v.re);
    let frame = owner
        .step_into(&mut state, input, &mut workspace, &mut || false)
        .unwrap();
    // Independent direct physical x/v generator, not the canonical pHS matrix.
    let expected_x = [
        8.956082781636439e-8,
        -8.108881188361974e-8,
        1.1681275232622287e-7,
    ];
    let expected_v = [
        0.0006192603832895811,
        -0.00042595093244568075,
        0.0007669034604454789,
    ];
    let expected_wall = [
        0.0009073019666416205,
        -0.00042595093244568075,
        0.0003852628729626785,
    ];
    let expected_traction = [-5355.860081893602, 1802.7728589334636, -3006.0875626920803];
    for axis in 0..3 {
        close(
            state.principal_displacement_m()[axis],
            expected_x[axis],
            5e-17,
        );
        close(
            state.principal_velocity_m_s()[axis],
            expected_v[axis],
            5e-14,
        );
        close(frame.support_velocity_m_s[axis], expected_wall[axis], 5e-14);
        close(
            frame.support_traction_pa[axis],
            expected_traction[axis],
            5e-7,
        );
    }
    close(
        frame.work.incoming_wave_energy_j_m2,
        0.0014304961136094675,
        3e-12,
    );
    close(
        frame.work.outgoing_wave_energy_j_m2,
        0.00025300083559744887,
        3e-12,
    );
    close(
        frame.work.support_dissipation_j_m2,
        0.0004934828864970197,
        3e-12,
    );
    close(
        frame.work.stored_energy_change_j_m2,
        0.0006840123915149989,
        3e-12,
    );
    assert!(frame.work.balance_residual_j_m2().abs() <= frame.work.balance_tolerance_j_m2);
}

#[test]
fn g3_native_normal_modes_recover_legacy_isotropic_fields_and_memory() {
    let material = solid(1.0);
    let characteristic = NormalElasticCharacteristic::new(&material, IDENTITY).unwrap();
    let support = support();
    let cartesian = [
        C64::from_re(0.0004),
        C64::from_re(-0.0002),
        C64::from_re(0.001),
    ];
    let input = modal(&characteristic, cartesian);
    let legacy = OrientedOscillatorSupportedSurface::prepare(
        &material,
        support,
        2000.0,
        ElasticMode::P,
        1.0,
    )
    .unwrap();
    let native =
        NormalOscillatorSupportedSurface::prepare(&characteristic, support, 2000.0).unwrap();
    let a = legacy
        .drive([cartesian[2], cartesian[0], cartesian[1]])
        .unwrap();
    let b = native.drive(input).unwrap();
    for axis in 0..3 {
        cc(
            a.reflection.boundary_velocity_m_s[axis],
            b.reflection.boundary_velocity_m_s[axis],
            5e-15,
        );
        cc(
            a.reflection.boundary_traction_pa[axis],
            b.reflection.boundary_traction_pa[axis],
            5e-8,
        );
    }
    close(a.mean_stored_energy_j_m2, b.mean_stored_energy_j_m2, 3e-13);
    let old_time = PreparedTransientElasticReflection::prepare(&material, support, 0.0002).unwrap();
    let new_time = PreparedTransientElasticReflection::prepare_normal_characteristic(
        &characteristic,
        support,
        0.0002,
    )
    .unwrap();
    let mut old_state = old_time.state(0.0, [1e-8; 3], [2e-5; 3]).unwrap();
    let mut new_state = old_state.clone();
    let mut old_work = old_time.workspace();
    let mut new_work = new_time.workspace();
    let a = old_time
        .step_into(
            &mut old_state,
            [0.001, 0.0004, -0.0002],
            &mut old_work,
            &mut || false,
        )
        .unwrap();
    let b = new_time
        .step_into(
            &mut new_state,
            input.map(|v| v.re),
            &mut new_work,
            &mut || false,
        )
        .unwrap();
    for axis in 0..3 {
        close(
            old_state.principal_velocity_m_s()[axis],
            new_state.principal_velocity_m_s()[axis],
            5e-14,
        );
        close(
            old_state.principal_displacement_m()[axis],
            new_state.principal_displacement_m()[axis],
            5e-17,
        );
        close(
            a.support_traction_pa[axis],
            b.support_traction_pa[axis],
            5e-7,
        );
    }
    close(
        a.work.outgoing_wave_energy_j_m2,
        b.work.outgoing_wave_energy_j_m2,
        3e-12,
    );
    close(
        a.work.support_dissipation_j_m2,
        b.work.support_dissipation_j_m2,
        3e-12,
    );
}

#[test]
fn g4_material_surface_orientation_binds_causal_state_transactionally() {
    let material = anisotropic();
    let characteristic = NormalElasticCharacteristic::new(&material, IDENTITY).unwrap();
    let other_frame = [
        [15.0 / 17.0, 0.0, 8.0 / 17.0],
        [0.0, 1.0, 0.0],
        [-8.0 / 17.0, 0.0, 15.0 / 17.0],
    ];
    let tilted = NormalElasticCharacteristic::new(&material, other_frame).unwrap();
    let owner = PreparedTransientElasticReflection::prepare_normal_characteristic(
        &characteristic,
        support(),
        0.0002,
    )
    .unwrap();
    let other = PreparedTransientElasticReflection::prepare_normal_characteristic(
        &tilted,
        support(),
        0.0002,
    )
    .unwrap();
    let mut state = owner.state(0.0, [1e-8; 3], [2e-5; 3]).unwrap();
    let before = state.clone();
    let mut workspace = owner.workspace();
    assert!(
        other
            .step_into(&mut state, [0.0; 3], &mut workspace, &mut || false)
            .is_err()
    );
    assert_eq!(state.time_s(), before.time_s());
    assert_eq!(
        state.principal_displacement_m(),
        before.principal_displacement_m()
    );
    assert_eq!(
        state.principal_velocity_m_s(),
        before.principal_velocity_m_s()
    );
    assert!(
        PreparedTransientElasticReflection::prepare_normal_characteristic_controlled(
            &characteristic,
            support(),
            0.0002,
            Default::default(),
            &mut || true,
        )
        .is_err()
    );
    assert!(
        NormalOscillatorSupportedSurface::prepare_controlled(
            &characteristic,
            support(),
            2000.0,
            &mut || true,
        )
        .is_err()
    );
    assert!(
        owner
            .step_into(&mut state, [0.0; 3], &mut workspace, &mut || false)
            .is_ok()
    );
}
