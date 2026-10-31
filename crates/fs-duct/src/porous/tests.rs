use super::*;
use crate::compose::{ComposeOptions, compose_chain, compose_scattering};
use crate::{ConeSlices, Duct, WallModel};
use fs_la::eigen_complex::lu_complex;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::PorousCard;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn material() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap()
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual {actual:?}, expected {expected:?}"
    );
}
fn phase(k: C64, length: f64) -> C64 {
    let phase = k.re * length;
    C64::new(phase.cos(), phase.sin()).scale((-k.im * length).exp())
}

// Independent four-interface-equation solve. Backward pressure is referenced
// at x=0 here, unlike the production stable endpoint bases. Only modest optical
// depths use this oracle: its growing exponential is unsuitable for production.
fn slab_reference(z: C64, k: C64, length: f64, reference: f64) -> (C64, C64) {
    let e = phase(k, length);
    let inverse = e.recip();
    let eta = C64::from_re(reference) / z;
    let zero = C64::ZERO;
    let one = C64::ONE;
    // Unknowns: reflected pressure, transmitted pressure, forward/backward
    // internal pressure. Incoming pressure is 1 Pa from the left, none right.
    let a = [
        one,
        zero,
        -one,
        -one,
        -one,
        zero,
        -eta,
        eta,
        zero,
        one,
        -e,
        -inverse,
        zero,
        one,
        -eta * e,
        eta * inverse,
    ];
    let mut rhs = [-one, -one, zero, zero];
    lu_complex(&a, 4).unwrap().solve(&mut rhs);
    (rhs[0], rhs[1])
}

#[test]
fn g1_slab_matches_pressure_and_superficial_velocity_continuity() {
    let gas = gas();
    let mut card = PorousCard::melamine_uf_nguyen2024();
    // Synthetic numerical fixture makes an erroneous second porosity factor
    // unmistakable. It is not a new measured material or source qualification.
    card.name = "synthetic-porosity-continuity";
    card.provenance = "synthetic algebra fixture, not measured material";
    card.porosity = 0.63;
    let material = PorousFluidLaw::try_new(card, 1.81e-5).unwrap();
    for hz in [200.0, 1000.0, 3500.0] {
        for area in [0.002, 0.008] {
            let omega = core::f64::consts::TAU * hz;
            let reference = gas.characteristic_impedance / area;
            let layer =
                PorousLayer::try_new(&material, &gas, omega, 0.04, area, reference).unwrap();
            assert!(std::ptr::eq(layer.material(), &material));
            assert_eq!(
                *layer.material_response().fluid(),
                gas.thermoviscous_state().unwrap()
            );
            assert_eq!(layer.area_m2(), area);
            assert_eq!(layer.thickness_m(), 0.04);
            let (specific, k) = layer.material_response().characteristic();
            let (r, t) =
                slab_reference(specific.conj().scale(1.0 / area), k.conj(), 0.04, reference);
            let power = layer
                .power_for_incident_pressure([C64::ONE, C64::ZERO])
                .unwrap();
            close(power.outgoing_pressure_pa[0], r, 3e-12);
            close(power.outgoing_pressure_pa[1], t, 3e-12);
            let absorption = power.absorbed_power_w / power.incident_power_w;
            assert!(absorption > 0.0 && absorption < 1.0);
            assert!((r.norm_sq() + t.norm_sq() + absorption - 1.0).abs() < 3e-12);
            assert!(
                t.norm_sq() > 1e-5,
                "this fixture must distinguish absorption from transmission"
            );
            assert!(((1.0 - r.norm_sq()) - absorption - t.norm_sq()).abs() < 3e-12);
        }
    }
}

#[test]
fn g1_coherent_power_matches_original_viscous_and_thermal_volume_loss() {
    let gas = gas();
    let material = material();
    let omega = core::f64::consts::TAU * 1100.0;
    let (area, length) = (0.003, 0.06);
    let layer = PorousLayer::try_new(
        &material,
        &gas,
        omega,
        length,
        area,
        gas.characteristic_impedance / area,
    )
    .unwrap();
    let sample = layer.material_response();
    let (z_plus, k_plus) = sample.characteristic();
    let (z, k) = (z_plus.conj(), k_plus.conj());
    let viscous = omega * sample.dynamic_density_kg_m3().conj().im;
    let thermal = omega * sample.dynamic_bulk_modulus_pa().conj().recip().im;
    assert!(viscous > 0.0 && thermal > 0.0);
    let drives = [
        [C64::ONE, C64::ZERO],
        [C64::ONE, C64::ONE],
        [C64::ONE, -C64::ONE],
        [C64::new(1.0, 0.3), C64::new(0.8, -0.2)],
    ];
    let mut powers = Vec::new();
    for drive in drives {
        let result = layer.power_for_incident_pressure(drive).unwrap();
        // Recover the independent physical field from the two face traces.
        let left = (result.pressure_pa[0]
            + z * result.inward_volume_flow_m3_s[0].scale(1.0 / area))
        .scale(0.5);
        let right = (result.pressure_pa[1]
            + z * result.inward_volume_flow_m3_s[1].scale(1.0 / area))
        .scale(0.5);
        let mut integral = 0.0;
        let n = 512;
        for j in 0..=n {
            let x = length * j as f64 / n as f64;
            let forward = left * phase(k, x);
            let backward = right * phase(k, length - x);
            let pressure = forward + backward;
            let superficial_velocity = (forward - backward) / z;
            let density = 0.5
                * area
                * (viscous * superficial_velocity.norm_sq() + thermal * pressure.norm_sq());
            integral += density
                * if j == 0 || j == n {
                    1.0
                } else if j % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
        }
        integral *= length / (3.0 * n as f64);
        assert!((integral - result.absorbed_power_w).abs() / integral < 2e-9);
        let boundary_power: f64 = result
            .pressure_pa
            .into_iter()
            .zip(result.inward_volume_flow_m3_s)
            .map(|(p, u)| 0.5 * (p * u.conj()).re)
            .sum();
        assert!((boundary_power - result.absorbed_power_w).abs() / integral < 2e-12);
        assert!(result.absorbed_power_w > 0.0);
        powers.push(result.absorbed_power_w);
    }
    assert!(
        (powers[1] - powers[2]).abs() > 0.01 * powers[1],
        "two equal-energy coherent drives must retain their different standing-wave absorption"
    );
}

fn loaded_slab(z: C64, k: C64, length: f64, load: C64) -> C64 {
    let e = phase(k, length);
    let cos = (e + e.recip()).scale(0.5);
    let sin = (e - e.recip()) / C64::new(0.0, 2.0);
    let minus_i = C64::new(0.0, -1.0);
    (load * cos + minus_i * z * sin) / (cos + minus_i * (load / z) * sin)
}

#[test]
fn g3_split_layers_rigid_backing_and_air_gap_use_original_composer() {
    let gas = gas();
    let material = material();
    let (omega, length, area) = (core::f64::consts::TAU * 900.0, 0.05, 0.004);
    let reference = gas.characteristic_impedance / area;
    let layer = PorousLayer::try_new(&material, &gas, omega, length, area, reference).unwrap();
    let half = PorousLayer::try_new(&material, &gas, omega, length / 2.0, area, reference).unwrap();
    let options = ComposeOptions::default();
    let (split, _) =
        compose_scattering(&[half.two_port(), half.two_port()], options, &mut || false).unwrap();
    for (actual, expected) in split
        .scattering
        .into_iter()
        .zip(layer.two_port().scattering)
    {
        close(actual, expected, 3e-12);
    }
    let closed = compose_chain(&[layer.two_port()], None, options).unwrap();
    close(
        closed.z_in.scale(area),
        layer
            .material_response()
            .surface_impedance_rigid(length)
            .unwrap()
            .conj(),
        3e-12,
    );
    assert_eq!(closed.u_out_over_u_in, C64::ZERO);
    // An actual inviscid gap uses the existing fluid duct chain, so the new
    // material layer can compose with original geometry and terminal loads.
    let gap_length = 0.04;
    let gap = Duct {
        segments: vec![Segment::Cylinder {
            radius: (area / core::f64::consts::PI).sqrt(),
            length: gap_length,
        }],
    };
    let gap = crate::fluid::duct_chain(
        &gap,
        &gas,
        omega,
        WallModel::Rigid,
        ConeSlices::Auto,
        options,
        || false,
    )
    .unwrap();
    let combined = compose_chain(&[layer.two_port(), gap.blocks[0]], None, options).unwrap();
    let gap_load = C64::new(
        0.0,
        reference / (omega * gap_length / gas.sound_speed).tan(),
    );
    let (z, k) = layer.material_response().characteristic();
    close(
        combined.z_in,
        loaded_slab(z.conj().scale(1.0 / area), k.conj(), length, gap_load),
        3e-12,
    );
    assert!((combined.z_in - closed.z_in).abs() > 0.1 * closed.z_in.abs());
    assert!(matches!(
        compose_scattering(&[layer.two_port()], options, &mut || true),
        Err(DuctError::Cancelled { completed_steps: 0 })
    ));
}

#[test]
fn g0_opaque_thin_zero_drive_and_source_refusals_remain_finite() {
    let gas = gas();
    let material = material();
    let (omega, area) = (core::f64::consts::TAU * 1000.0, 0.002);
    let reference = gas.characteristic_impedance / area;
    let opaque = PorousLayer::try_new(&material, &gas, omega, 1000.0, area, reference).unwrap();
    let z = opaque
        .material_response()
        .characteristic()
        .0
        .conj()
        .scale(1.0 / area);
    let half_space = (z - C64::from_re(reference)) / (z + C64::from_re(reference));
    close(opaque.two_port().scattering[0], half_space, 3e-12);
    assert_eq!(opaque.two_port().scattering[1], C64::ZERO);
    assert_eq!(opaque.two_port().scattering[2], C64::ZERO);
    let power = opaque
        .power_for_incident_pressure([C64::ONE, C64::ZERO])
        .unwrap();
    assert!(power.absorbed_power_w > 0.0);
    let zero = opaque.power_for_incident_pressure([C64::ZERO; 2]).unwrap();
    assert_eq!(zero.absorbed_power_w, 0.0);
    assert_eq!(zero.pressure_pa, [C64::ZERO; 2]);
    let thin = PorousLayer::try_new(&material, &gas, omega, 1e-10, area, reference).unwrap();
    close(thin.two_port().scattering[0], C64::ZERO, 1e-8);
    close(thin.two_port().scattering[2], C64::ONE, 1e-8);
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(PorousLayer::try_new(&material, &gas, omega, invalid, area, reference).is_err());
        assert!(PorousLayer::try_new(&material, &gas, omega, 0.05, invalid, reference).is_err());
        assert!(PorousLayer::try_new(&material, &gas, omega, 0.05, area, invalid).is_err());
    }
    assert!(PorousLayer::try_new(&material, &gas, 1.0, 0.05, area, reference).is_err());
    let mut stale = gas;
    stale.pressure *= 2.0;
    assert!(PorousLayer::try_new(&material, &stale, omega, 0.05, area, reference).is_err());
    assert!(
        opaque
            .power_for_incident_pressure([C64::new(f64::NAN, 0.0), C64::ZERO])
            .is_err()
    );
    assert!(
        opaque
            .power_for_incident_pressure([C64::from_re(f64::MAX); 2])
            .is_err()
    );
}
