//! G1/G3/G4: original anisotropic coating fields and actual rear pin storage.

use super::{IDENTITY, anisotropic, cc, close, modal, solid, support};
use fs_blake3::ContentHash;
use fs_couple::elastic_reflection::{
    ElasticSupportError, NormalOscillatorBackedElasticCoating,
    NormalOscillatorSupportedSurface, OrientedOscillatorSupport, OscillatorBackedElasticCoating,
};
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::interface::ElasticMode;
use fs_material::elastic::acoustic::interface::normal::NormalElasticCharacteristic;
use fs_math::c64::C64;

fn film(limit: f64) -> ElasticAcousticState {
    // Independent normal eigenvalues rho*c^2 for rho=1800 and
    // c=[1100,2100,3400], with noncommuting mixed polarizations.
    let mut c = [[0.0; 6]; 6];
    for (i, value) in [
        16e9,
        14e9,
        18_052_082_840.23669,
        9e9,
        12_725_463_656.25594,
        13_018_370_663.270693,
    ]
    .into_iter()
    .enumerate()
    {
        c[i][i] = value;
    }
    for (i, j, value) in [
        (2, 4, 6_773_508_737.0809145),
        (2, 5, 6_450_960_701.981824),
        (4, 5, -3_000_510_803.565776),
    ] {
        c[i][j] = value;
        c[j][i] = value;
    }
    ElasticAcousticState::try_new_mandel(c, 1800.0, ContentHash([57; 32]), limit).unwrap()
}

fn input() -> [C64; 3] {
    [
        C64::from_re(0.001),
        C64::new(-0.0002, 0.0001),
        C64::new(0.0003, -0.0004),
    ]
}

#[test]
fn g3_absent_anisotropic_coating_recovers_original_direct_support_exactly() {
    let incident = anisotropic();
    let coating = film(0.002);
    let a = NormalElasticCharacteristic::new(&incident, IDENTITY).unwrap();
    let b = NormalElasticCharacteristic::new(&coating, IDENTITY).unwrap();
    let frequency = 2000.0 / core::f64::consts::TAU;
    let omega = core::f64::consts::TAU * frequency;
    let direct = NormalOscillatorSupportedSurface::prepare(&a, support(), omega).unwrap();
    let coated = NormalOscillatorBackedElasticCoating::prepare(
        &a, &b, 0.0, support(), frequency,
    )
    .unwrap();
    let incoming = modal(&a, input());
    let direct = direct.drive(incoming).unwrap();
    let response = coated.drive(incoming).unwrap();
    assert!(core::ptr::eq(
        coated.coating().incident_characteristic().state(),
        &incident
    ));
    assert!(core::ptr::eq(
        coated.coating().coating_characteristic().state(),
        &coating
    ));
    assert_eq!(coated.support(), &support());
    assert_eq!(coated.omega_rad_s(), omega);
    assert_eq!(coated.coating().frequency_hz(), frequency);
    assert!(response.coating.coating_faces.is_none());
    assert!(response.coating.field_at(0.0).is_err());
    assert_eq!(
        response.coating.reflection.reflected_velocity_m_s,
        direct.reflection.reflected_velocity_m_s
    );
    assert_eq!(response.principal_velocity_m_s, direct.principal_velocity_m_s);
    assert_eq!(response.support_displacement_m, direct.support_displacement_m);
    assert_eq!(
        response.mean_support_stored_energy_j_m2,
        direct.mean_stored_energy_j_m2
    );
}

#[test]
fn g1_anisotropic_coating_and_rear_pin_storage_match_independent_cartesian_transfer() {
    let incident = anisotropic();
    let coating = film(0.002);
    let a = NormalElasticCharacteristic::new(&incident, IDENTITY).unwrap();
    let b = NormalElasticCharacteristic::new(&coating, IDENTITY).unwrap();
    let owner = NormalOscillatorBackedElasticCoating::prepare(
        &a,
        &b,
        0.43,
        support(),
        2000.0 / core::f64::consts::TAU,
    )
    .unwrap();
    let response = owner.drive_cartesian(input()).unwrap();
    // Independent binary64 Cartesian continuum transfer matrix and 3x3
    // physical boundary solve; it does not reuse the native modal layer solve.
    let reflected = [
        C64::new(0.0003366177109050661, -0.0001014314239409185),
        C64::new(-0.0002912990675367962, 0.00004852294981711704),
        C64::new(0.00011238885506990482, -0.0002607171153104454),
    ];
    let rear = [
        C64::new(0.0009950075348079768, -0.00005844025006261862),
        C64::new(-0.0005424828164728077, -0.000317229437654918),
        C64::new(0.0005506830519272021, -0.00033126713749752676),
    ];
    let traction = [
        C64::new(-3269.806797788582, -2513.663153192891),
        C64::new(1944.6778870733413, 409.2054964919463),
        C64::new(-4387.742268698551, 990.0892283745795),
    ];
    let principal = [
        C64::new(0.0006188028592177667, 0.000104325491120055),
        C64::new(-0.0005424828164728077, -0.000317229437654918),
        C64::new(0.000954135650433638, -0.00031979582723316765),
    ];
    let kinetic = [
        0.00009845019666838096,
        0.00009873053057078137,
        0.00025316105263604186,
    ];
    let elastic = [
        0.0002953505900051429,
        0.00004936526528539068,
        0.0005063221052720837,
    ];
    for i in 0..3 {
        cc(response.coating.reflection.reflected_velocity_m_s[i], reflected[i], 2e-13);
        cc(response.coating.support_velocity_m_s[i], rear[i], 2e-13);
        cc(response.coating.support_traction_pa[i], traction[i], 3e-6);
        cc(response.principal_velocity_m_s[i], principal[i], 2e-13);
        cc(
            response.principal_displacement_m[i] * C64::new(0.0, -owner.omega_rad_s()),
            principal[i],
            2e-13,
        );
        close(response.mean_kinetic_energy_j_m2[i], kinetic[i], 5e-13);
        close(response.mean_elastic_energy_j_m2[i], elastic[i], 5e-13);
    }
    close(
        response.mean_support_stored_energy_j_m2,
        0.0013013797404378214,
        8e-13,
    );
    close(response.coating.reflection.total_incoming_power_w_m2, 4.428288284023669, 3e-9);
    close(response.coating.reflection.total_reflected_power_w_m2, 0.9104937875650404, 3e-9);
    close(response.coating.reflection.absorbed_normal_power_w_m2, 3.51779449645863, 3e-9);
    for depth in [0.0, 0.07, 0.21, 0.43] {
        close(
            response.coating.field_at(depth).unwrap().normal_flux_w_m2,
            3.51779449645863,
            3e-9,
        );
    }
    // Rear motion, not incident-face motion, owns support energy.
    let front = response.coating.reflection.boundary_velocity_m_s;
    assert!(
        front.into_iter().zip(rear).any(|(front, rear)| (front - rear).abs() > 1e-5)
    );
}

#[test]
fn g3_native_anisotropic_api_recovers_existing_isotropic_coating_fields_and_storage() {
    let incident = solid(1.0);
    let coating = solid(0.7);
    let a = NormalElasticCharacteristic::new(&incident, IDENTITY).unwrap();
    let b = NormalElasticCharacteristic::new(&coating, IDENTITY).unwrap();
    let frequency = 2000.0 / core::f64::consts::TAU;
    let normal = NormalOscillatorBackedElasticCoating::prepare(
        &a, &b, 0.43, support(), frequency,
    )
    .unwrap();
    let legacy = OscillatorBackedElasticCoating::prepare(
        &incident,
        &coating,
        0.43,
        support(),
        frequency,
        ElasticMode::P,
        1.0,
    )
    .unwrap();
    let cartesian = input();
    let a = normal.drive_cartesian(cartesian).unwrap();
    let b = legacy.drive([cartesian[2], cartesian[0], cartesian[1]]).unwrap();
    for i in 0..3 {
        cc(a.coating.support_velocity_m_s[i], b.coating.support_velocity_m_s[i], 2e-13);
        cc(a.coating.support_traction_pa[i], b.coating.support_traction_pa[i], 3e-6);
        cc(a.principal_displacement_m[i], b.principal_displacement_m[i], 2e-16);
    }
    for depth in [0.0, 0.15, 0.43] {
        let af = a.coating.field_at(depth).unwrap();
        let bf = b.coating.field_at(depth).unwrap();
        for i in 0..3 {
            cc(af.velocity_m_s[i], bf.velocity_m_s[i], 2e-13);
            cc(af.traction_pa[i], bf.traction_pa[i], 3e-6);
        }
    }
    close(a.mean_support_stored_energy_j_m2, b.mean_support_stored_energy_j_m2, 8e-13);
    close(
        a.coating.reflection.absorbed_normal_power_w_m2,
        b.coating.reflection.absorbed_normal_power_w_m2,
        3e-9,
    );
}

#[test]
fn g1_lossless_coating_and_reactive_support_store_energy_without_heat() {
    let incident = anisotropic();
    let coating = film(0.002);
    let a = NormalElasticCharacteristic::new(&incident, IDENTITY).unwrap();
    let b = NormalElasticCharacteristic::new(&coating, IDENTITY).unwrap();
    let original = support();
    let pins = (*original.principal_supports()).map(|mut pin| {
        pin.resistance = 0.0;
        pin
    });
    let reactive = OrientedOscillatorSupport::new(pins, original.principal_to_surface());
    let owner = NormalOscillatorBackedElasticCoating::prepare(
        &a,
        &b,
        0.43,
        reactive,
        2000.0 / core::f64::consts::TAU,
    )
    .unwrap();
    let response = owner.drive_cartesian(input()).unwrap();
    assert_eq!(response.coating.reflection.absorbed_normal_power_w_m2, 0.0);
    close(
        response.coating.reflection.total_incoming_power_w_m2,
        response.coating.reflection.total_reflected_power_w_m2,
        3e-9,
    );
    assert!(response.mean_support_stored_energy_j_m2 > 0.0);
}

#[test]
fn g4_original_coating_domain_and_cancelled_preparation_do_not_publish_fields() {
    let incident = anisotropic();
    let coating = film(1e-12);
    let a = NormalElasticCharacteristic::new(&incident, IDENTITY).unwrap();
    let b = NormalElasticCharacteristic::new(&coating, IDENTITY).unwrap();
    let frequency = 2000.0 / core::f64::consts::TAU;
    let owner = NormalOscillatorBackedElasticCoating::prepare(
        &a, &b, 0.43, support(), frequency,
    )
    .unwrap();
    assert!(owner.drive_cartesian(input()).is_err());
    assert!(owner.drive([C64::ZERO; 3]).is_ok());
    assert!(matches!(
        NormalOscillatorBackedElasticCoating::prepare_controlled(
            &a, &b, 0.43, support(), frequency, &mut || true,
        ),
        Err(ElasticSupportError::Cancelled)
    ));
    assert!(
        NormalOscillatorBackedElasticCoating::prepare(
            &a, &b, -0.43, support(), frequency,
        )
        .is_err()
    );
    assert!(
        NormalOscillatorBackedElasticCoating::prepare(
            &a, &b, 0.43, support(), f64::NAN,
        )
        .is_err()
    );
}
