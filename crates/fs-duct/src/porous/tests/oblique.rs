use super::*;
use crate::compose::ComposeOptions;

#[test]
fn g1_oblique_slab_matches_independent_interface_continuity_and_exact_normal_limit() {
    let gas = gas();
    let material = material();
    let (omega, area, length) = (core::f64::consts::TAU * 1200.0, 0.004, 0.035);
    let reference = gas.characteristic_impedance / area;
    let normal = PorousLayer::try_new(&material, &gas, omega, length, area, reference).unwrap();
    let zero =
        PorousLayer::try_new_oblique(&material, &gas, omega, length, area, reference, 0.0).unwrap();
    assert_eq!(normal.two_port().scattering, zero.two_port().scattering);
    let drive = [C64::new(1.0, 0.2), C64::new(0.3, -0.6)];
    let normal_fields = normal.prepare_excitation(drive).unwrap();
    let zero_fields = zero.prepare_excitation(drive).unwrap();
    assert_eq!(
        normal_fields.at_depth(0.01).unwrap(),
        zero_fields.at_depth(0.01).unwrap()
    );
    assert_eq!(
        normal_fields.dissipation_between(0.0, length).unwrap(),
        zero_fields.dissipation_between(0.0, length).unwrap()
    );
    for cosine in [0.95_f64, 0.6, 0.08] {
        let kt = omega / gas.sound_speed * (1.0 - cosine * cosine).sqrt();
        let reference = reference / cosine;
        let layer =
            PorousLayer::try_new_oblique(&material, &gas, omega, length, area, reference, kt)
                .unwrap();
        let response = layer.material_response();
        let rho = response.dynamic_density_kg_m3().conj();
        let bulk = response.dynamic_bulk_modulus_pa().conj();
        let q = ((rho / bulk).scale(omega * omega) - C64::from_re(kt * kt)).sqrt();
        let z = rho.scale(omega) / q;
        close(layer.normal_wavenumber_rad_m(), q, 3e-14);
        close(layer.normal_specific_impedance_pa_s_m(), z, 3e-14);
        let (r, t) = slab_reference(z.scale(1.0 / area), q, length, reference);
        close(layer.two_port().scattering[0], r, 2e-12);
        close(layer.two_port().scattering[2], t, 2e-12);
        assert_ne!(normal.two_port().scattering, layer.two_port().scattering);
    }
}

#[test]
fn g1_transverse_velocity_adds_positive_viscous_work_and_closes_normal_flux() {
    let gas = gas();
    let material = material();
    let (omega, area, length, cosine) = (6000.0, 0.004, 0.06, 0.35_f64);
    let kt = omega / gas.sound_speed * (1.0 - cosine * cosine).sqrt();
    let reference = gas.characteristic_impedance / (area * cosine);
    let layer =
        PorousLayer::try_new_oblique(&material, &gas, omega, length, area, reference, kt).unwrap();
    let excitation = layer
        .prepare_excitation([C64::new(1.0, 0.3), C64::new(-0.2, 0.6)])
        .unwrap();
    let response = layer.material_response();
    let rho = response.dynamic_density_kg_m3().conj();
    let bulk = response.dynamic_bulk_modulus_pa().conj();
    let q = ((rho / bulk).scale(omega * omega) - C64::from_re(kt * kt)).sqrt();
    let z = rho.scale(omega) / q;
    let boundary = excitation.boundary_power();
    let p0 = boundary.pressure_pa[0];
    let v0 = boundary.inward_volume_flow_m3_s[0].scale(1.0 / area);
    let forward = (p0 + z * v0).scale(0.5);
    let backward = (p0 - z * v0).scale(0.5);
    let mut integrals = [0.0; 3];
    let n = 2048;
    for i in 0..=n {
        let x = length * i as f64 / n as f64;
        let f = forward * phase(q, x);
        let b = backward * phase(-q, x);
        let p = f + b;
        let vn = (f - b) / z;
        let vt = p * (C64::from_re(kt / omega) / rho);
        let viscous_n = 0.5 * omega * rho.im * vn.norm_sq();
        let viscous_t = 0.5 * omega * rho.im * vt.norm_sq();
        let thermal = 0.5 * omega * bulk.recip().im * p.norm_sq();
        let field = excitation.at_depth(x).unwrap();
        close(field.pressure_pa, p, 5e-12);
        close(field.superficial_velocity_m_s, vn, 5e-12);
        close(field.tangential_superficial_velocity_m_s, vt, 5e-12);
        assert!((field.viscous_dissipation_w_m3 - viscous_n - viscous_t).abs() < 5e-13);
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, value) in integrals.iter_mut().zip([viscous_n, viscous_t, thermal]) {
            *sum += weight * value;
        }
    }
    let integrated = integrals.map(|x| x * area * length / (3.0 * n as f64));
    let work = excitation.dissipation_between(0.0, length).unwrap();
    assert!(
        integrated[1] > 0.01 * work.absorbed_power_w,
        "the transverse term must be material"
    );
    assert!((work.viscous_power_w - integrated[0] - integrated[1]).abs() < 2e-14);
    assert!((work.thermal_power_w - integrated[2]).abs() < 2e-14);
    assert!((work.absorbed_power_w - boundary.absorbed_power_w).abs() < 2e-14);
}

#[test]
fn g0_oblique_mode_support_evanescence_and_opaque_limits_remain_bounded() {
    let gas = gas();
    let material = material();
    let (omega, area, reference) = (6000.0, 0.004, 100_000.0);
    let normal = PorousLayer::try_new(&material, &gas, omega, 0.04, area, reference).unwrap();
    let k = normal.material_response().characteristic().1.conj();
    for kt in [k.re * (1.0 - 1e-10), k.re * (1.0 + 1e-10), 3.0 * k.abs()] {
        let layer = PorousLayer::try_new_oblique(&material, &gas, omega, 0.04, area, reference, kt)
            .unwrap();
        assert!(layer.normal_wavenumber_rad_m().im > 0.0);
        let excitation = layer
            .prepare_excitation([C64::ONE, C64::new(0.2, 0.4)])
            .unwrap();
        let work = excitation.dissipation_between(0.0, 0.04).unwrap();
        assert!(work.absorbed_power_w > 0.0);
        assert!(
            (work.absorbed_power_w - excitation.boundary_power().absorbed_power_w).abs() < 2e-14
        );
        let mismatch = [normal.clone(), layer];
        assert!(PorousLayerStack::try_new(&mismatch, ComposeOptions::default(), || false).is_err());
    }
    let opaque =
        PorousLayer::try_new_oblique(&material, &gas, omega, 1000.0, area, reference, k.re)
            .unwrap();
    let dark = opaque.prepare_excitation([C64::ZERO; 2]).unwrap();
    assert_eq!(
        dark.dissipation_between(0.0, 1000.0)
            .unwrap()
            .absorbed_power_w,
        0.0
    );
    let fields = opaque.prepare_excitation([C64::ONE, C64::ZERO]).unwrap();
    assert_eq!(
        fields
            .at_depth(500.0)
            .unwrap()
            .tangential_superficial_velocity_m_s,
        C64::ZERO
    );
    assert_eq!(opaque.two_port().scattering[2], C64::ZERO);
    for kt in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PorousLayer::try_new_oblique(&material, &gas, omega, 0.04, area, reference, kt)
                .is_err()
        );
    }
}
