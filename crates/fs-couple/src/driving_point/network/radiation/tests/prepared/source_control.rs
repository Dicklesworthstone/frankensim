use super::*;
use fs_bem::observation_transfers::{
    BemObservation, HarmonicControlLimits, HarmonicSourceMapLimits, HarmonicVolumeQuadrature,
    SourceId,
};

const MAP: HarmonicSourceMapLimits = HarmonicSourceMapLimits {
    maximum_sources: 32,
    maximum_products: 100_000,
    maximum_storage_bytes: 100_000,
};
const CONTROL: HarmonicControlLimits = HarmonicControlLimits {
    maximum_quadrature_points: 16,
    maximum_products: 100_000,
    maximum_workspace_bytes: 100_000,
};
const VOLUMES: [f64; 2] = [0.03, 0.08];

fn quadrature() -> [HarmonicVolumeQuadrature; 2] {
    core::array::from_fn(|receiver| HarmonicVolumeQuadrature {
        receiver,
        volume_m3: VOLUMES[receiver],
    })
}
fn energy(medium: Medium, fields: &[HarmonicField], targets: Option<&[HarmonicField]>) -> f64 {
    fields
        .iter()
        .enumerate()
        .map(|(r, field)| {
            let target = targets.map_or(HarmonicField::default(), |t| t[r]);
            let pressure = field.pressure_pa - target.pressure_pa;
            let velocity = (0..3)
                .map(|axis| {
                    (field.particle_velocity_m_s[axis] - target.particle_velocity_m_s[axis])
                        .norm_sq()
                })
                .sum::<f64>();
            VOLUMES[r]
                * (pressure.norm_sq() / (medium.density * medium.sound_speed.powi(2))
                    + medium.density * velocity)
                * 0.25
        })
        .sum()
}
fn fresh_fields(owner: &PreparedRadiatingNetwork<'_>, source: &[C64]) -> Vec<HarmonicField> {
    let gate = CancelGate::new_clock_free();
    let mut response = owner.response_buffer();
    owner.drive_into(source, &mut response, &gate).unwrap();
    direct_fields(
        owner.model(),
        owner.omega(),
        response.outward_flows_m3_s().unwrap(),
    )
}
fn scalar_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-23),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn loaded_mixed_si_generator_gradients_match_fresh_physical_fields() {
    let model = two_source_model();
    let omega = core::f64::consts::TAU * 1200.0;
    let sources = [
        PortSource {
            label: "inlet".into(),
            excitation: PortExcitation::VolumeFlow(C64::new(0.8e-6, -0.3e-6)),
        },
        PortSource {
            label: "aux".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: C64::new(0.6, 0.2),
                impedance_pa_s_m3: C64::new(2e6, -1e6),
            },
        },
    ];
    let gate = CancelGate::new_clock_free();
    let owner = model
        .prepare_sources(omega, LossModel::Bessel, &sources, None, &[], &gate)
        .unwrap();
    let bridge = owner.source_control(MAP, &gate).unwrap();
    assert!(std::ptr::eq(bridge.owner(), &owner));
    assert_eq!(bridge.source_coordinates(), owner.source_coordinates());
    assert_eq!(bridge.default_excitations(), owner.default_excitations());
    assert_eq!(
        bridge.observation().sources(),
        &[
            SourceId::OpenCircuitPressure("aux".into()),
            SourceId::PrescribedVolumeFlow("inlet".into()),
        ]
    );
    assert!(owner.boundary().z[1].abs() > 0.0);
    let x = bridge.default_excitations().to_vec();
    let direct = fresh_fields(&owner, &x);
    let mut observed = vec![HarmonicField::default(); direct.len()];
    bridge
        .observation()
        .evaluate_coherent(&x, &mut observed, || false)
        .unwrap();
    for (a, b) in observed.iter().zip(&direct) {
        close(a.pressure_pa, b.pressure_pa);
        for axis in 0..3 {
            close(a.particle_velocity_m_s[axis], b.particle_velocity_m_s[axis]);
        }
    }
    let target = [HarmonicField {
        pressure_pa: C64::new(1.0e-4, -2.0e-4),
        particle_velocity_m_s: [C64::new(2.0e-7, -1.0e-7); 3],
    }; 2];
    let volume = quadrature();
    let medium = Medium::from(&model.gas());
    let mut control = bridge
        .energy_control(&volume, Some(&target), CONTROL, &gate)
        .unwrap();
    let result = control.value_gradient(&x, || false).unwrap();
    scalar_close(
        result.mean_residual_energy_j,
        energy(medium, &direct, Some(&target)),
        2e-8,
    );
    let gradient = result.source_gradient.to_vec();
    let output_address = result.source_gradient.as_ptr();
    // Independent fresh BEM boundary fields, with SI steps for Pa and m³/s.
    for source in 0..2 {
        let step = if source == 0 { 1.0e-5 } else { 1.0e-11 };
        for imaginary in [false, true] {
            let delta = if imaginary {
                C64::new(0.0, step)
            } else {
                C64::from_re(step)
            };
            let (mut plus, mut minus) = (x.clone(), x.clone());
            plus[source] = plus[source] + delta;
            minus[source] = minus[source] - delta;
            let finite_difference = (energy(medium, &fresh_fields(&owner, &plus), Some(&target))
                - energy(medium, &fresh_fields(&owner, &minus), Some(&target)))
                / (2.0 * step);
            scalar_close(
                if imaginary {
                    gradient[source].im
                } else {
                    gradient[source].re
                },
                finite_difference,
                3e-5,
            );
        }
    }
    let direction = [C64::new(0.3, -0.6), C64::new(0.2e-6, 0.5e-6)];
    let hv = control.hessian_vector(&direction, || false).unwrap();
    assert_eq!(hv.as_ptr(), output_address);
    let directional = direction
        .iter()
        .zip(hv)
        .map(|(&a, &b)| (a.conj() * b).re)
        .sum::<f64>();
    scalar_close(
        directional,
        2.0 * energy(medium, &fresh_fields(&owner, &direction), None),
        3e-8,
    );
    // No extra radiation work is charged by preparing/evaluating controls.
    let mut after = owner.response_buffer();
    owner.drive_into(&x, &mut after, &gate).unwrap();
    scalar_close(
        after.network().unwrap().fields.external_load_power_w,
        after.exterior_surface_power_w().unwrap(),
        3e-9,
    );
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        owner.source_control(MAP, &stopped),
        Err(ExteriorNetworkError::Cancelled)
    ));
    assert!(matches!(
        control.value_gradient(&x, || true),
        Err(fs_bem::observation_transfers::TransferError::Cancelled)
    ));
    assert_eq!(
        control
            .value_gradient(&x, || false)
            .unwrap()
            .source_gradient
            .as_ptr(),
        output_address
    );
}

#[test]
fn mutual_generator_impedance_remains_upstream_of_energy_controls_and_trace_admission() {
    let model = two_source_model();
    let omega = core::f64::consts::TAU * 1200.0;
    let coupled = bank(omega);
    let mut diagonal = bank(omega);
    diagonal.impedance.z[1] = C64::ZERO;
    diagonal.impedance.z[2] = C64::ZERO;
    let gate = CancelGate::new_clock_free();
    let owner = model
        .prepare_sources(omega, LossModel::Lossless, &[], Some(&coupled), &[], &gate)
        .unwrap();
    let uncoupled = model
        .prepare_sources(omega, LossModel::Lossless, &[], Some(&diagonal), &[], &gate)
        .unwrap();
    let bridge = owner.source_control(MAP, &gate).unwrap();
    let other = uncoupled.source_control(MAP, &gate).unwrap();
    assert_ne!(
        bridge.observation().producer_hash(),
        other.observation().producer_hash()
    );
    assert!(
        bridge
            .observation()
            .sources()
            .iter()
            .all(|s| matches!(s, SourceId::OpenCircuitPressure(_)))
    );
    let volume = quadrature();
    let mut control = bridge
        .energy_control(&volume, None, CONTROL, &gate)
        .unwrap();
    let mut other_control = other.energy_control(&volume, None, CONTROL, &gate).unwrap();
    let x = owner.default_excitations();
    let value = control
        .value_gradient(x, || false)
        .unwrap()
        .mean_residual_energy_j;
    let other_value = other_control
        .value_gradient(x, || false)
        .unwrap()
        .mean_residual_energy_j;
    scalar_close(
        value,
        energy(Medium::from(&model.gas()), &fresh_fields(&owner, x), None),
        2e-8,
    );
    assert!((value - other_value).abs() > 1.0e-4 * value.abs());
    // Upstream generators cannot be smuggled in as surface velocity traces.
    for kind in [
        SourceId::OpenCircuitPressure("generator".into()),
        SourceId::PrescribedVolumeFlow("generator".into()),
    ] {
        assert!(
            BemObservation::new(
                model.surface(),
                Medium::from(&model.gas()),
                &[kind],
                &[vec![1.0; model.surface().areas().len()]],
                bridge.observation().receivers(),
                0.001
            )
            .is_err()
        );
    }
}
