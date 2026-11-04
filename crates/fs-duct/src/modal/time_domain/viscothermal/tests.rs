//! G1/G2/G3: physical wall work, mode traces, and independent Bessel propagation.
use super::*;
use fs_material::gas::GasSpec;
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel};

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn budget() -> ModalViscothermalBudget {
    ModalViscothermalBudget {
        maximum_states: 124,
        maximum_payload_bytes: 32 * 1024 * 1024,
        maximum_preparation_work: 200_000_000,
    }
}
fn resolution(
    gas: &GasState,
    radius: f64,
    shear: f64,
    poles: usize,
) -> ModalViscothermalResolution {
    let omega = shear * shear * gas.dynamic_viscosity / (gas.density * radius * radius);
    ModalViscothermalResolution {
        poles,
        band_rad_s: [0.99 * omega, 1.01 * omega],
        maximum_diffusive_relative_error: 0.5,
        maximum_boundary_layer_ratio: 0.1,
    }
}
fn make(elements: usize, modes: usize, end: ModalEnd, poles: usize, shear: f64) -> ModalGalerkin {
    let gas = air();
    let profile = RadiusProfile::exponential(0.002, 0.0, 0.006).unwrap();
    assemble_lossless(&profile, &gas, ModeCount(modes), elements, end, 124, || {
        false
    })
    .unwrap()
    .with_viscothermal(resolution(&gas, 0.002, shear, poles), budget(), || false)
    .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

/// Populate ORIGINAL histories with their steady phasors for independently
/// supplied physical q and P. No harmonic equation or R is used here.
fn harmonic_state(
    source: &ModalGalerkinSource,
    q: &[C64],
    pressure: &[C64],
    omega: f64,
) -> Vec<C64> {
    let n = source.coordinates();
    let wall = source.viscothermal.as_ref().unwrap();
    let mut state = vec![C64::ZERO; source.state_dimension()];
    state[..n].copy_from_slice(q);
    for i in 0..n {
        state[n + i] = (0..n).fold(C64::ZERO, |sum, j| {
            sum + pressure[j].scale(source.mass[i * n + j])
        });
    }
    let mut cursor = 2 * n;
    for (channels, field, terms) in [
        (&wall.velocity, q, wall.model.series_terms()),
        (&wall.pressure, pressure, wall.model.shunt_terms()),
    ] {
        for term in terms {
            let coupling = det::sqrt(term.gain_per_s * term.rate_per_s);
            for row in channels {
                let value = row
                    .iter()
                    .zip(field)
                    .fold(C64::ZERO, |sum, (&a, &x)| sum + x.scale(a));
                state[cursor] = value.scale(coupling) / C64::new(term.rate_per_s, -omega);
                cursor += 1;
            }
        }
    }
    state
}

#[test]
fn higher_mode_wall_work_uses_actual_signed_pressure_and_tangential_velocity() {
    let gas = air();
    let gamma = Basis::new(2).gammas[1];
    // Both modes propagate at this actual frequency; beta_1/beta_0=0.6.
    let omega = gas.sound_speed * gamma / (0.8 * 0.002);
    let shear = 0.002 * det::sqrt(omega * gas.density / gas.dynamic_viscosity);
    let image = make(1, 2, ModalEnd::Ports, 38, shear);
    let source = image.source();
    let n = source.coordinates();
    let gas = source.gas();
    let radius = source.profile.inlet_radius();
    let length = source.profile.length();
    let omega = source.viscothermal_model().unwrap().band_rad_s()[0] / 0.99;
    // Independent exact Bessel surface impedance/admittance. The dimensions
    // here are wall area, local velocity and pressure, not canonical damping.
    let rv = radius * det::sqrt(omega * gas.density / gas.dynamic_viscosity);
    let s = C64::new(0.0, omega);
    let shear_impedance =
        (s / fs_phs::zk_one_minus_f(rv).unwrap() - s).scale(gas.density * radius / 2.0);
    let thermal_admittance = (s * fs_phs::zwikker_kosten_f(rv * det::sqrt(gas.prandtl)).unwrap())
        .scale((gas.gamma - 1.0) * radius / (2.0 * gas.density * gas.sound_speed.powi(2)));
    let wall_area = core::f64::consts::TAU * radius * length;
    let plane_gradient = 3.0 / gas.sound_speed;
    let mut heats = Vec::new();
    for mode in 0..2 {
        let gradient = plane_gradient * if mode == 0 { 1.0 } else { 0.6 };
        let mut q = vec![C64::ZERO; n];
        let mut pressure = vec![C64::ZERO; n];
        for (node, &x) in source.nodes.iter().enumerate() {
            q[2 * node + mode] = C64::from_re(gradient * x);
            pressure[2 * node + mode] = C64::from_re(3.0);
        }
        let state = harmonic_state(source, &q, &pressure, omega);
        let heat = source.viscothermal_mean_heat_w(&state).unwrap();
        close(
            heat.viscous_w,
            0.5 * wall_area * shear_impedance.re * (gradient / gas.density).powi(2),
            0.001,
        );
        close(
            heat.thermal_w,
            0.5 * wall_area * thermal_admittance.re * 9.0,
            0.001,
        );
        heats.push(heat);
    }
    // At a common frequency a propagating mode with beta_1/beta_0=0.6
    // has 0.36 of the shear work for equal pressure, while thermal work agrees.
    close(heats[1].viscous_w / heats[0].viscous_w, 0.36, 2e-12);
    close(heats[1].thermal_w / heats[0].thermal_w, 1.0, 2e-12);
    // psi_0(R)=+1, psi_1(R)=-1: equal modal amplitudes cancel at the
    // wall while retaining nonzero bulk energy. Independent scalar mode loss
    // cannot reproduce this interference or its sign-reversed enhancement.
    for sign in [1.0, -1.0] {
        let mut q = vec![C64::ZERO; n];
        let mut p = vec![C64::ZERO; n];
        for (node, &x) in source.nodes.iter().enumerate() {
            q[2 * node] = C64::from_re(plane_gradient * x);
            q[2 * node + 1] = C64::from_re(sign * plane_gradient * x);
            p[2 * node] = C64::from_re(3.0);
            p[2 * node + 1] = C64::from_re(sign * 3.0);
        }
        let state = harmonic_state(source, &q, &p, omega);
        let heat = source.viscothermal_mean_heat_w(&state).unwrap();
        if sign == 1.0 {
            assert!(heat.viscous_w < 1e-25 * heats[0].viscous_w);
            assert!(heat.thermal_w < 1e-25 * heats[0].thermal_w);
        } else {
            close(heat.viscous_w, 4.0 * heats[0].viscous_w, 2e-12);
            close(heat.thermal_w, 4.0 * heats[0].thermal_w, 2e-12);
        }
    }
}

#[test]
fn modal_wall_power_approaches_independent_flat_boundary_layer_law() {
    let mut errors = Vec::new();
    for (shear, poles) in [(100.0, 16), (200.0, 24), (400.0, 38)] {
        let image = make(1, 2, ModalEnd::Ports, poles, shear);
        let source = image.source();
        let gas = source.gas();
        let n = source.coordinates();
        let omega = source.viscothermal_model().unwrap().band_rad_s()[0] / 0.99;
        let mut q = vec![C64::ZERO; n];
        let mut p = vec![C64::ZERO; n];
        for (node, &x) in source.nodes.iter().enumerate() {
            q[2 * node + 1] = C64::from_re(x);
            p[2 * node + 1] = C64::ONE;
        }
        let heat = source
            .viscothermal_mean_heat_w(&harmonic_state(source, &q, &p, omega))
            .unwrap();
        let area = core::f64::consts::TAU * 0.002 * source.profile.length();
        // Independently solved oscillating half-space Stokes/Fourier layers:
        // Re Z_shear=sqrt(rho*mu*omega/2), Re Y_T=(gamma-1)/(rho*c²)*sqrt(alpha*omega/2).
        let viscous = 0.5 * area * det::sqrt(gas.density * gas.dynamic_viscosity * omega / 2.0)
            / gas.density.powi(2);
        let alpha = gas.dynamic_viscosity / (gas.density * gas.prandtl);
        let thermal = 0.5 * area * (gas.gamma - 1.0) / (gas.density * gas.sound_speed.powi(2))
            * det::sqrt(alpha * omega / 2.0);
        errors.push([
            (heat.viscous_w / viscous - 1.0).abs(),
            (heat.thermal_w / thermal - 1.0).abs(),
        ]);
    }
    for mechanism in 0..2 {
        assert!(errors[0][mechanism] > 1e-3, "{errors:?}");
        assert!(
            errors[1][mechanism] < 0.65 * errors[0][mechanism],
            "{errors:?}"
        );
        assert!(
            errors[2][mechanism] < 0.65 * errors[1][mechanism],
            "{errors:?}"
        );
        assert!(errors[2][mechanism] < 0.012, "{errors:?}");
    }
}

#[test]
fn original_diffusion_poles_refine_both_physical_wall_losses() {
    let mut errors = Vec::new();
    for poles in [8, 12, 16] {
        let image = make(1, 2, ModalEnd::Ports, poles, 100.0);
        let source = image.source();
        let model = source.viscothermal_model().unwrap();
        errors.push([model.max_complex_error(), model.max_real_loss_error()]);
        assert_eq!(source.viscothermal_state_range().unwrap().len(), 3 * poles);
        assert!(image.system.resistance_witness().is_some());
    }
    for quantity in 0..2 {
        assert!(errors[0][quantity] > 0.1, "{errors:?}");
        assert!(
            errors[1][quantity] < 0.1 * errors[0][quantity],
            "{errors:?}"
        );
        assert!(
            errors[2][quantity] < 0.05 * errors[1][quantity],
            "{errors:?}"
        );
    }
}

fn harmonic_budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 124,
        max_dense_bytes: 32 * 1024 * 1024,
        max_dense_work: 8_000_000,
        relative_residual_tolerance: 2e-10,
    }
}

#[test]
fn plane_complete_wall_state_converges_to_independent_bessel_duct() {
    let gas = air();
    let r = 0.002;
    let resolution = resolution(&gas, r, 100.0, 16);
    let omega = resolution.band_rad_s[0] / 0.99;
    let length = 0.7 * gas.sound_speed / omega;
    let profile = RadiusProfile::exponential(r, 0.0, length).unwrap();
    let duct = crate::Duct {
        segments: vec![crate::Segment::Cylinder { radius: r, length }],
    };
    let zref = gas.density * gas.sound_speed / (core::f64::consts::PI * r * r);
    for (end, termination) in [
        (ModalEnd::Rigid, crate::Termination::Closed),
        (ModalEnd::PressureRelease, crate::Termination::IdealOpen),
    ] {
        let exact = crate::input_impedance(&duct, &gas, omega, LossModel::Bessel, termination)
            .unwrap()
            .impedance;
        let reference = C64::from_re(zref);
        let exact_reflection = (exact - reference) / (exact + reference);
        let mut errors = Vec::new();
        for elements in [1, 2, 3] {
            let image =
                assemble_lossless(&profile, &gas, ModeCount(1), elements, end, 124, || false)
                    .unwrap()
                    .with_viscothermal(resolution, budget(), || false)
                    .unwrap();
            let model =
                HarmonicModel::from_ode(image.system(), harmonic_budget(), || false).unwrap();
            let response = model
                .factor(omega, || false)
                .unwrap()
                .response(&[C64::ONE], || false)
                .unwrap();
            let z = response.output[0];
            errors.push(((z - reference) / (z + reference) - exact_reflection).abs());
            let heat = image
                .source()
                .viscothermal_mean_heat_w(&response.state)
                .unwrap();
            close(
                heat.viscous_w + heat.thermal_w,
                response.dissipated_power_w,
                2e-8,
            );
            close(response.supplied_power_w, response.dissipated_power_w, 2e-8);
            assert!(heat.viscous_w > 0.0 && heat.thermal_w > 0.0);
        }
        assert!(errors[0] > 1e-4, "{errors:?}");
        assert!(
            errors[1] < 0.35 * errors[0] && errors[2] < 0.55 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 0.008, "{errors:?}");
    }
}

#[test]
fn complete_two_face_response_is_reciprocal_and_heat_closes_port_work() {
    let image = make(2, 2, ModalEnd::Ports, 16, 100.0);
    let source = image.source();
    let omega = source.viscothermal_model().unwrap().band_rad_s()[0] / 0.99;
    let model = HarmonicModel::from_ode(image.system(), harmonic_budget(), || false).unwrap();
    let factor = model.factor(omega, || false).unwrap();
    let mut z = [C64::ZERO; 16];
    for column in 0..4 {
        let mut drive = [C64::ZERO; 4];
        drive[column] = C64::ONE;
        let response = factor.response(&drive, || false).unwrap();
        for row in 0..4 {
            z[row * 4 + column] = response.output[row];
        }
        let heat = source.viscothermal_mean_heat_w(&response.state).unwrap();
        close(
            heat.viscous_w + heat.thermal_w,
            response.dissipated_power_w,
            2e-8,
        );
        close(response.supplied_power_w, response.dissipated_power_w, 2e-8);
    }
    let scale = source.inlet_reference_impedance();
    for i in 0..4 {
        for j in 0..4 {
            assert!((z[i * 4 + j] - z[j * 4 + i]).abs() < 2e-9 * scale);
        }
    }
    assert!(
        z[1].abs() > 1e-5 * scale,
        "physical wall cross-mode coupling must be live"
    );
}

#[test]
fn exact_gauge_and_original_coordinates_survive_both_end_conditions() {
    for end in [ModalEnd::Rigid, ModalEnd::PressureRelease, ModalEnd::Ports] {
        let image = make(2, 2, end, 16, 100.0);
        let source = image.source();
        assert_eq!(source.state_dimension(), image.system.state_dim());
        let mut state = vec![0.0; source.state_dimension()];
        let n = source.coordinates();
        for i in 0..n {
            state[i] = det::sin(i as f64);
            state[n + i] = 1e-10 * det::cos(i as f64);
        }
        for (i, value) in state.iter_mut().enumerate().skip(2 * n) {
            *value = 1e-4 * det::sin(i as f64);
        }
        let heat = source.viscothermal_heat_w(&state).unwrap();
        let effort = image.system.effort(&state);
        let r = image.system.structure().1;
        let power: f64 = (0..state.len())
            .map(|i| {
                (0..state.len())
                    .map(|j| effort[i] * r[i * state.len() + j] * effort[j])
                    .sum::<f64>()
            })
            .sum();
        close(heat.viscous_w + heat.thermal_w, power, 2e-10);
        if end != ModalEnd::PressureRelease {
            for node in 0..source.nodes.len() {
                state[node * 2] += 7.0;
            }
            let shifted = source.viscothermal_heat_w(&state).unwrap();
            close(shifted.viscous_w, heat.viscous_w, 2e-12);
            assert_eq!(shifted.thermal_w, heat.thermal_w);
        }
    }
    // A perturbed physical stiffness must not be repaired as if it were a gauge.
    let mut k = vec![1.01, -1.0, -1.0, 1.0];
    assert!(lift_velocity(&k, 1, ModalEnd::Ports, &[vec![-1.0, 1.0]], &mut || false).is_err());
    k[0] = 1.0;
    assert!(lift_velocity(&k, 1, ModalEnd::Ports, &[vec![-1.0, 1.0]], &mut || false).is_ok());
}

#[test]
fn original_shape_capacity_thin_layer_and_cancellation_refusals_remain_explicit() {
    let gas = air();
    let shape = RadiusProfile::exponential(0.002, 0.0, 0.006).unwrap();
    let resolution = resolution(&gas, 0.002, 100.0, 16);
    let bare = || {
        assemble_lossless(&shape, &gas, ModeCount(2), 2, ModalEnd::Ports, 124, || {
            false
        })
        .unwrap()
    };
    let mut small = budget();
    small.maximum_states = 91;
    assert!(
        bare()
            .with_viscothermal(resolution, small, || false)
            .is_err()
    );
    let mut narrow = resolution;
    narrow.band_rad_s = resolution.band_rad_s.map(|f| 0.1 * f);
    assert!(
        bare()
            .with_viscothermal(narrow, budget(), || false)
            .is_err()
    );
    let tapered = RadiusProfile::exponential(0.002, 0.1, 0.006).unwrap();
    assert!(
        assemble_lossless(
            &tapered,
            &gas,
            ModeCount(2),
            2,
            ModalEnd::Ports,
            124,
            || false
        )
        .unwrap()
        .with_viscothermal(resolution, budget(), || false)
        .is_err()
    );
    assert!(
        make(2, 2, ModalEnd::Ports, 16, 100.0)
            .with_viscothermal(resolution, budget(), || false)
            .is_err()
    );
    let medium = fs_material::relaxation::RelaxingMedium::new(
        gas,
        vec![(
            fs_thermochem::SpeciesId::new("authored-empty-law").unwrap(),
            1.0,
        )],
    )
    .unwrap();
    let molecular = fs_material::relaxation::RelaxationAbsorption::new(vec![])
        .unwrap()
        .evaluate(&medium)
        .unwrap();
    let molecular_budget = ModalMolecularBudget {
        maximum_states: 124,
        maximum_payload_bytes: 32 * 1024 * 1024,
        maximum_preparation_work: 200_000_000,
    };
    assert!(
        bare()
            .with_molecular_relaxation(molecular.clone(), molecular_budget, || false)
            .unwrap()
            .with_viscothermal(resolution, budget(), || false)
            .is_err()
    );
    assert!(
        make(2, 2, ModalEnd::Ports, 16, 100.0)
            .with_molecular_relaxation(molecular, molecular_budget, || false)
            .is_err()
    );
    let mut total = 0;
    bare()
        .with_viscothermal(resolution, budget(), || {
            total += 1;
            false
        })
        .unwrap();
    for limit in [1, 7, total] {
        let mut polls = 0;
        let result = bare().with_viscothermal(resolution, budget(), || {
            polls += 1;
            polls == limit
        });
        assert!(matches!(
            result,
            Err(ModalTimeError::Duct(DuctError::Cancelled { .. }))
        ));
        assert_eq!(polls, limit);
    }
    assert_eq!(shape.constant_radius_m(), Some(0.002));
    assert_eq!(
        RadiusProfile::power(0.002, 1.0, 0.0, 0.006)
            .unwrap()
            .constant_radius_m(),
        Some(0.002)
    );
    let stations = vec![(0.0, 0.002), (0.003, 0.002), (0.006, 0.002)];
    assert_eq!(
        RadiusProfile::spline(stations.clone(), None)
            .unwrap()
            .constant_radius_m(),
        Some(0.002)
    );
    assert_eq!(
        RadiusProfile::spline(stations, Some([0.01, -0.01]))
            .unwrap()
            .constant_radius_m(),
        None
    );
}
