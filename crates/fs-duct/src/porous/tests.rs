use super::*;
use crate::compose::{ComposeOptions, compose_chain, compose_scattering};
use crate::{ConeSlices, Duct, WallModel};
use fs_la::eigen_complex::lu_complex;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::PorousCard;

mod oblique;

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
        let excitation = layer.prepare_excitation(drive).unwrap();
        assert!(std::ptr::eq(excitation.layer(), &layer));
        let result = excitation.boundary_power();
        // Recover the independent physical field from the two face traces.
        let left = (result.pressure_pa[0]
            + z * result.inward_volume_flow_m3_s[0].scale(1.0 / area))
        .scale(0.5);
        let right = (result.pressure_pa[1]
            + z * result.inward_volume_flow_m3_s[1].scale(1.0 / area))
        .scale(0.5);
        let mut viscous_integral = 0.0;
        let mut thermal_integral = 0.0;
        let n = 512;
        for j in 0..=n {
            let x = length * j as f64 / n as f64;
            let forward = left * phase(k, x);
            let backward = right * phase(k, length - x);
            let pressure = forward + backward;
            let superficial_velocity = (forward - backward) / z;
            let field = excitation.at_depth(x).unwrap();
            close(field.pressure_pa, pressure, 3e-12);
            close(field.superficial_velocity_m_s, superficial_velocity, 3e-12);
            close(
                field.volume_flow_m3_s,
                superficial_velocity.scale(area),
                3e-12,
            );
            let viscous_density = 0.5 * viscous * superficial_velocity.norm_sq();
            let thermal_density = 0.5 * thermal * pressure.norm_sq();
            assert!((field.viscous_dissipation_w_m3 - viscous_density).abs() < 2e-14);
            assert!((field.thermal_dissipation_w_m3 - thermal_density).abs() < 2e-14);
            let weight = if j == 0 || j == n {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            viscous_integral += area * viscous_density * weight;
            thermal_integral += area * thermal_density * weight;
        }
        viscous_integral *= length / (3.0 * n as f64);
        thermal_integral *= length / (3.0 * n as f64);
        let integral = viscous_integral + thermal_integral;
        let exact = excitation.dissipation_between(0.0, length).unwrap();
        assert!((exact.viscous_power_w / viscous_integral - 1.0).abs() < 2e-9);
        assert!((exact.thermal_power_w / thermal_integral - 1.0).abs() < 2e-9);
        assert!((exact.absorbed_power_w / result.absorbed_power_w - 1.0).abs() < 2e-12);
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

#[test]
fn g3_thermal_cells_preserve_power_partition_and_coherent_node_heating() {
    let gas = gas();
    let material = material();
    let (area, length) = (0.003, 0.06);
    for hz in [200.0, 1100.0, 3500.0] {
        let layer = PorousLayer::try_new(
            &material,
            &gas,
            core::f64::consts::TAU * hz,
            length,
            area,
            gas.characteristic_impedance / area,
        )
        .unwrap();
        for drive in [
            [C64::ONE, C64::ONE],
            [C64::ONE, -C64::ONE],
            [C64::new(0.7, 0.3), C64::new(-0.2, 0.8)],
        ] {
            let excitation = layer.prepare_excitation(drive).unwrap();
            let whole = excitation.dissipation_between(0.0, length).unwrap();
            let mut viscous = 0.0;
            let mut thermal = 0.0;
            for j in 0..37 {
                let start = length * j as f64 / 37.0;
                let end = length * (j + 1) as f64 / 37.0;
                let cell = excitation.dissipation_between(start, end).unwrap();
                assert!(cell.viscous_power_w >= 0.0 && cell.thermal_power_w >= 0.0);
                viscous += cell.viscous_power_w;
                thermal += cell.thermal_power_w;
                let a = excitation.at_depth(start).unwrap();
                let b = excitation.at_depth(end).unwrap();
                let flux_difference = 0.5
                    * ((a.pressure_pa * a.volume_flow_m3_s.conj()).re
                        - (b.pressure_pa * b.volume_flow_m3_s.conj()).re);
                assert!((cell.absorbed_power_w - flux_difference).abs() < 2e-17);
            }
            assert!((viscous / whole.viscous_power_w - 1.0).abs() < 3e-13);
            assert!((thermal / whole.thermal_power_w - 1.0).abs() < 3e-13);
            let left = excitation.at_depth(0.0).unwrap();
            let right = excitation.at_depth(length).unwrap();
            close(
                left.pressure_pa,
                excitation.boundary_power().pressure_pa[0],
                3e-12,
            );
            close(
                right.pressure_pa,
                excitation.boundary_power().pressure_pa[1],
                3e-12,
            );
            close(
                right.volume_flow_m3_s,
                -excitation.boundary_power().inward_volume_flow_m3_s[1],
                3e-12,
            );
        }
        // A symmetric drive has a velocity node at the midpoint. Viscous
        // heating in a tiny centered cell is still positive (order width^3).
        let excitation = layer.prepare_excitation([C64::ONE; 2]).unwrap();
        let width = 1e-9;
        let center = 0.5 * length;
        let cell = excitation
            .dissipation_between(center - 0.5 * width, center + 0.5 * width)
            .unwrap();
        assert!(cell.viscous_power_w > 0.0);
        let pressure = excitation.at_depth(center).unwrap().pressure_pa;
        let response = layer.material_response();
        let velocity_gradient = pressure * C64::new(0.0, response.angular_frequency_rad_s())
            / response.dynamic_bulk_modulus_pa().conj();
        let leading = 0.5
            * response.angular_frequency_rad_s()
            * -response.dynamic_density_kg_m3().im
            * area
            * velocity_gradient.norm_sq()
            * width.powi(3)
            / 12.0;
        assert!((cell.viscous_power_w / leading - 1.0).abs() < 2e-8);
    }
}

#[test]
fn g0_pore_heating_respects_zero_thermal_coupling_and_extensive_scaling() {
    use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
    let gas = gas();
    let original = gas.thermoviscous_state().unwrap();
    // Manufactured general-fluid Cp=Cv, retaining original c and rho. This
    // admits no thermal acoustic coupling and is not a new gas EOS.
    let fluid = ThermoviscousFluidState::from_heat_capacity_ratio(
        original.acoustic(),
        original.transport(),
        1.0,
    )
    .unwrap();
    let material = material();
    let mut powers = Vec::new();
    for area in [0.002, 0.006] {
        let layer = PorousLayer::try_new(
            &material,
            &fluid,
            6000.0,
            0.08,
            area,
            gas.characteristic_impedance / area,
        )
        .unwrap();
        let excitation = layer.prepare_excitation([C64::ONE, C64::ZERO]).unwrap();
        assert_eq!(
            excitation.at_depth(0.03).unwrap().thermal_dissipation_w_m3,
            0.0
        );
        let heat = excitation.dissipation_between(0.0, 0.08).unwrap();
        assert_eq!(heat.thermal_power_w, 0.0);
        assert_eq!(heat.absorbed_power_w, heat.viscous_power_w);
        assert!(
            (heat.absorbed_power_w / excitation.boundary_power().absorbed_power_w - 1.0).abs()
                < 3e-12
        );
        powers.push(heat.absorbed_power_w);
    }
    assert!((powers[1] / powers[0] - 3.0).abs() < 3e-12);
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
    let excitation = opaque.prepare_excitation([C64::ONE, C64::ZERO]).unwrap();
    let heat = excitation.dissipation_between(0.0, 1000.0).unwrap();
    assert!((heat.absorbed_power_w / power.absorbed_power_w - 1.0).abs() < 3e-12);
    let middle = excitation.at_depth(500.0).unwrap();
    assert_eq!(middle.pressure_pa, C64::ZERO);
    assert_eq!(middle.viscous_dissipation_w_m3, 0.0);
    assert_eq!(middle.thermal_dissipation_w_m3, 0.0);
    for invalid in [-0.1, 1000.1, f64::NAN, f64::INFINITY] {
        assert!(excitation.at_depth(invalid).is_err());
        assert!(excitation.dissipation_between(0.0, invalid).is_err());
        assert!(excitation.dissipation_between(invalid, 1000.0).is_err());
    }
    assert!(excitation.dissipation_between(1.0, 0.0).is_err());
    assert_eq!(
        excitation
            .dissipation_between(0.2, 0.2)
            .unwrap()
            .absorbed_power_w,
        0.0
    );
    let zero = opaque.power_for_incident_pressure([C64::ZERO; 2]).unwrap();
    assert_eq!(zero.absorbed_power_w, 0.0);
    assert_eq!(zero.pressure_pa, [C64::ZERO; 2]);
    assert_eq!(
        opaque
            .prepare_excitation([C64::ZERO; 2])
            .unwrap()
            .dissipation_between(0.0, 1000.0)
            .unwrap()
            .absorbed_power_w,
        0.0,
    );
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

fn molecular_source(active: bool, time_scale: f64) -> RelaxationState {
    use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
    // Manufactured reservoir for a constitutive test, not measured air data.
    let text = if active {
        format!(
            "declared 1 synthetic 1 vibration 1 900 1 power {} 293.15 0 200 500",
            101_325.0 * 0.001 * time_scale,
        )
    } else {
        "declared 1 synthetic 0".into()
    };
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let medium = RelaxingMedium::new(gas(), vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}

#[test]
fn molecular_porous_slab_keeps_original_dispersion_and_oblique_interfaces() {
    let source = molecular_source(true, 1.0);
    let material = material();
    let (area, length) = (0.004, 0.04);
    let mut maximum_contrast = 0.0_f64;
    for frequency in [200.0, 1100.0, 3500.0] {
        let omega = core::f64::consts::TAU * frequency;
        for sine in [0.0, 0.75] {
            let kt = (omega / source.gas().sound_speed) * sine;
            let reference =
                source.gas().characteristic_impedance / (area * (1.0_f64 - sine * sine).sqrt());
            let layer = PorousLayer::try_new_relaxing_oblique(
                &material, &source, omega, length, area, reference, kt,
            )
            .unwrap();
            assert!(std::ptr::eq(layer.material(), &material));
            assert!(std::ptr::eq(layer.relaxation_state().unwrap(), &source));
            let original = material.effective_relaxing(&source, omega).unwrap();
            assert_eq!(
                layer.material_response().dynamic_bulk_modulus_pa(),
                original.dynamic_bulk_modulus_pa()
            );
            assert_eq!(
                layer.material_response().dynamic_density_kg_m3(),
                original.dynamic_density_kg_m3()
            );
            let rho = original.dynamic_density_kg_m3().conj();
            let bulk = original.dynamic_bulk_modulus_pa().conj();
            let q = ((rho / bulk).scale(omega * omega) - C64::from_re(kt * kt)).sqrt();
            let z = rho.scale(omega) / q;
            close(layer.normal_wavenumber_rad_m(), q, 3e-13);
            close(layer.normal_specific_impedance_pa_s_m(), z, 3e-13);
            let (r, t) = slab_reference(z.scale(1.0 / area), q, length, reference);
            close(layer.two_port().scattering[0], r, 3e-12);
            close(layer.two_port().scattering[2], t, 3e-12);
            let plain = PorousLayer::try_new_oblique(
                &material,
                source.gas(),
                omega,
                length,
                area,
                reference,
                kt,
            )
            .unwrap();
            maximum_contrast = maximum_contrast.max((t - plain.two_port().scattering[2]).abs());
            assert!(plain.relaxation_state().is_none());
        }
    }
    assert!(
        maximum_contrast > 1e-5,
        "original molecular storage must reach finite transmission"
    );
    for frequency in [99.0, 4001.0] {
        assert!(
            PorousLayer::try_new_relaxing(
                &material,
                &source,
                core::f64::consts::TAU * frequency,
                length,
                area,
                100_000.0,
            )
            .is_err()
        );
    }
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PorousLayer::try_new_relaxing(&material, &source, 6000.0, invalid, area, 100_000.0)
                .is_err()
        );
        assert!(
            PorousLayer::try_new_relaxing(&material, &source, 6000.0, length, invalid, 100_000.0)
                .is_err()
        );
        assert!(
            PorousLayer::try_new_relaxing(&material, &source, 6000.0, length, area, invalid)
                .is_err()
        );
    }
    assert!(
        PorousLayer::try_new_relaxing_oblique(
            &material, &source, 6000.0, length, area, 100_000.0, -1.0,
        )
        .is_err()
    );
}

#[test]
fn molecular_porous_coherent_heat_matches_independent_oblique_volume_work() {
    let source = molecular_source(true, 1.0);
    let material = material();
    let (area, length) = (0.004, 0.05);
    let omega = core::f64::consts::TAU * 900.0;
    let kt = 0.6 * omega / source.gas().sound_speed;
    let reference = source.gas().characteristic_impedance / (0.8 * area);
    let layer = PorousLayer::try_new_relaxing_oblique(
        &material, &source, omega, length, area, reference, kt,
    )
    .unwrap();
    let rho = layer.material_response().dynamic_density_kg_m3().conj();
    let bulk = layer.material_response().dynamic_bulk_modulus_pa().conj();
    let q = ((rho / bulk).scale(omega * omega) - C64::from_re(kt * kt)).sqrt();
    let z = rho.scale(omega) / q;
    let mut absorbed = Vec::new();
    for incoming in [
        [C64::ONE; 2],
        [C64::ONE, -C64::ONE],
        [C64::new(0.7, 0.3), C64::new(-0.2, 0.8)],
    ] {
        let excitation = layer.prepare_excitation(incoming).unwrap();
        let boundary = excitation.boundary_power();
        let left = (boundary.pressure_pa[0]
            + z * boundary.inward_volume_flow_m3_s[0].scale(1.0 / area))
        .scale(0.5);
        let right = (boundary.pressure_pa[1]
            + z * boundary.inward_volume_flow_m3_s[1].scale(1.0 / area))
        .scale(0.5);
        let mut integral = [0.0; 2];
        let n = 512;
        for i in 0..=n {
            let x = length * i as f64 / n as f64;
            let forward = left * phase(q, x);
            let backward = right * phase(q, length - x);
            let pressure = forward + backward;
            let normal_velocity = (forward - backward) / z;
            let tangential_velocity = pressure / rho * C64::from_re(kt / omega);
            let field = excitation.at_depth(x).unwrap();
            close(field.pressure_pa, pressure, 3e-12);
            close(field.superficial_velocity_m_s, normal_velocity, 3e-12);
            close(
                field.tangential_superficial_velocity_m_s,
                tangential_velocity,
                3e-12,
            );
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral[0] += weight
                * 0.5
                * omega
                * rho.im
                * (normal_velocity.norm_sq() + tangential_velocity.norm_sq());
            integral[1] += weight * 0.5 * omega * bulk.recip().im * pressure.norm_sq();
        }
        let integral = integral.map(|value| value * area * length / (3.0 * n as f64));
        let native = excitation.dissipation_between(0.0, length).unwrap();
        assert!(native.viscous_power_w > 0.0 && native.thermal_power_w > 0.0);
        assert!((native.viscous_power_w / integral[0] - 1.0).abs() < 2e-9);
        assert!((native.thermal_power_w / integral[1] - 1.0).abs() < 2e-9);
        assert!((native.absorbed_power_w / boundary.absorbed_power_w - 1.0).abs() < 2e-11);
        let first = excitation.dissipation_between(0.0, 0.017).unwrap();
        let second = excitation.dissipation_between(0.017, length).unwrap();
        assert!(
            ((first.absorbed_power_w + second.absorbed_power_w) / native.absorbed_power_w - 1.0)
                .abs()
                < 3e-13
        );
        absorbed.push(native.absorbed_power_w);
    }
    assert!((absorbed[0] - absorbed[1]).abs() > 0.01 * absorbed[0]);
}

#[test]
fn molecular_porous_stacks_preserve_empty_sources_subdivision_and_rear_load_work() {
    let empty = molecular_source(false, 1.0);
    let active = molecular_source(true, 1.0);
    let changed = molecular_source(true, 4.0);
    let material = material();
    let (area, length) = (0.004, 0.06);
    let omega = core::f64::consts::TAU * 1100.0;
    let reference = empty.gas().characteristic_impedance / area;
    for kt in [0.0, 0.6 * omega / empty.gas().sound_speed] {
        let plain = PorousLayer::try_new_oblique(
            &material,
            empty.gas(),
            omega,
            length,
            area,
            reference,
            kt,
        )
        .unwrap();
        let retained = PorousLayer::try_new_relaxing_oblique(
            &material, &empty, omega, length, area, reference, kt,
        )
        .unwrap();
        assert!(std::ptr::eq(retained.relaxation_state().unwrap(), &empty));
        assert_eq!(retained.two_port().scattering, plain.two_port().scattering);
        let drive = [C64::new(0.7, -0.2), C64::new(-0.1, 0.5)];
        let a = plain.prepare_excitation(drive).unwrap();
        let b = retained.prepare_excitation(drive).unwrap();
        assert_eq!(a.at_depth(0.023).unwrap(), b.at_depth(0.023).unwrap());
        assert_eq!(
            a.dissipation_between(0.005, 0.04).unwrap(),
            b.dissipation_between(0.005, 0.04).unwrap()
        );
    }
    let make = |state, depth| {
        PorousLayer::try_new_relaxing(&material, state, omega, depth, area, reference).unwrap()
    };
    let one = [make(&active, length)];
    let split = [make(&active, 0.017), make(&active, length - 0.017)];
    let one_stack = PorousLayerStack::try_new(&one, ComposeOptions::default(), || false).unwrap();
    let split_stack =
        PorousLayerStack::try_new(&split, ComposeOptions::default(), || false).unwrap();
    for (a, b) in one_stack
        .two_port()
        .scattering
        .into_iter()
        .zip(split_stack.two_port().scattering)
    {
        close(a, b, 3e-12);
    }
    let changed_layer = make(&changed, length);
    assert!(
        (changed_layer.two_port().scattering[2] - one[0].two_port().scattering[2]).abs() > 1e-6
    );
    let incident = C64::new(0.7, 0.2);
    let load = C64::new(0.8 * reference, 0.3 * reference);
    let termination = PassivePorousTermination::try_new_impedance(load).unwrap();
    let a = one_stack
        .prepare_terminated_excitation(incident, termination, || false)
        .unwrap();
    let b = split_stack
        .prepare_terminated_excitation(incident, termination, || false)
        .unwrap();
    let (z, k) = one[0].material_response().characteristic();
    let expected_impedance = loaded_slab(z.conj().scale(1.0 / area), k.conj(), length, load);
    let expected_reflection = (expected_impedance - C64::from_re(reference))
        / (expected_impedance + C64::from_re(reference));
    close(
        a.boundary_power().outgoing_pressure_pa[0],
        expected_reflection * incident,
        3e-12,
    );
    close(
        b.boundary_power().outgoing_pressure_pa[0],
        expected_reflection * incident,
        3e-12,
    );
    for excitation in [&a, &b] {
        let power = excitation.source_power();
        assert!(power.load_absorbed_power_w > 0.0 && power.porous_absorbed_power_w > 0.0);
        assert!(power.residual_w.abs() < 3e-13 * power.incident_power_w);
        let rear = excitation.at_depth(length).unwrap().field;
        let load_work = 0.5 * load.re * rear.volume_flow_m3_s.norm_sq();
        assert!((power.load_absorbed_power_w / load_work - 1.0).abs() < 3e-12);
    }
    assert!(
        (a.dissipation().absorbed_power_w / b.dissipation().absorbed_power_w - 1.0).abs() < 3e-12
    );
    close(
        a.at_depth(0.031).unwrap().field.pressure_pa,
        b.at_depth(0.031).unwrap().field.pressure_pa,
        3e-12,
    );
    assert!(matches!(
        split_stack.prepare_terminated_excitation(incident, termination, || true),
        Err(DuctError::Cancelled { .. })
    ));
}
