use super::*;
use crate::{
    helmholtz::{Formulation, Medium, bc::Admittance, laplace::complex::ComplexOptions},
    observation_transfers::{HarmonicField, LocalBoundary, LocalWall, ReceiverSpec},
    panel3d::SpherePanels,
};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};
use fs_phs::WallPin;

fn material() -> RelaxationState {
    // Manufactured positive caloric reservoir; not an experimental air fit.
    let text = "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500";
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}

fn observer<'a>(
    surface: &'a SpherePanels,
    source: RelaxationState,
    ids: &[SourceId],
    traces: &[Vec<f64>],
) -> MolecularBemObservation<'a> {
    let p = surface.areas().len();
    MolecularBemObservation::new(
        surface,
        source,
        ids,
        traces,
        &[ReceiverSpec {
            id: "physical microphone".into(),
            position_m: [0.0, 0.0, 1.0],
        }],
        0.01,
        [1.0, 3_000.0],
        ComplexOptions {
            max_panels: p,
            max_bytes: 32 * p * p,
            max_work: 512 * (p * p) as u64,
        },
    )
    .unwrap()
}

fn limits(model: &MolecularBemObservation<'_>) -> MolecularResponseLimits {
    let s = model.boundary.sources.len();
    let p = model.boundary.surface.areas().len();
    let r = model.boundary.receivers.len();
    MolecularResponseLimits {
        maximum_retained_entries: 3 * s * s + 2 * s * p + 8 * s * r,
        maximum_panel_products: (3 * s * s * p) as u64,
    }
}

fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * scale.max(f64::MIN_POSITIVE),
        "{a} versus {b}"
    );
}

#[test]
fn loaded_spherical_oscillator_and_receiver_converge_to_independent_solution() {
    // G1: physical mechanical impedance plus the analytic pulsating-sphere
    // fluid impedance. The same computed velocity drives the receiver and
    // reaction; no independently adjusted radiation damping is introduced.
    let source = material();
    let omega = 600.0;
    let radius = 0.1;
    let gas = *source.gas();
    let k = source.relaxing_wavenumber(omega).unwrap();
    let area = 4.0 * core::f64::consts::PI * radius * radius;
    let specific = C64::new(0.0, gas.density * omega) / C64::new(-1.0 / radius - k.im, k.re);
    let mechanical = C64::new(0.2, 0.01 * (400.0_f64.powi(2) / omega - omega));
    let force = C64::new(0.3, 0.1);
    let expected_velocity = force / (mechanical + specific.scale(area));
    let phase = k.scale(1.0 - radius);
    let propagation = C64::new(phase.re.cos(), phase.re.sin()).scale(radius * (-phase.im).exp());
    let expected_pressure = specific * expected_velocity * propagation;
    let mut errors = Vec::new();
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(radius, level).unwrap();
        let model = observer(
            &surface,
            source.clone(),
            &[SourceId::StructuralMode("radial displacement".into())],
            &[vec![1.0; surface.areas().len()]],
        );
        let response = model
            .prepare_response(omega, limits(&model), || false)
            .unwrap();
        let velocity = force / (mechanical + response.impedance()[0]);
        let mut field = [HarmonicField::default()];
        response
            .harmonic()
            .evaluate_coherent(&[velocity], &mut field, || false)
            .unwrap();
        let power = response.power_w(&[velocity], || false).unwrap();
        let input = 0.5 * (force * velocity.conj()).re;
        let mechanical_heat = 0.5 * mechanical.re * velocity.norm_sq();
        close(input, mechanical_heat + power.driving_power_w, input.abs());
        assert_eq!(power.absorbed_wall_power_w, 0.0);
        assert_eq!(response.wall_work_matrix(), &[C64::ZERO]);
        close(power.exterior_power_w, power.driving_power_w, input.abs());
        assert!(power.exterior_power_w > 0.0);
        errors.push(
            ((velocity - expected_velocity).abs() / expected_velocity.abs())
                .max((field[0].pressure_pa - expected_pressure).abs() / expected_pressure.abs()),
        );
    }
    assert!(
        errors[2] < 0.12 && errors[2] < 0.65 * errors[0],
        "{errors:?}"
    );
}

#[test]
fn mixed_source_work_keeps_wall_absorption_and_coherent_mutual_loading() {
    // G0/G3: a mode and a unit-volume-flow opening drive the lower hemisphere;
    // the upper hemisphere is a passive physical mass/spring/dashpot wall.
    // Compare the cached quadratic matrices against an original-owner solve
    // of the complete complex drive, including receiver p/gradient.
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let source = material();
    let medium = Medium::from(source.gas());
    let mask: Vec<_> = surface.centroids().iter().map(|x| x[2] <= 0.0).collect();
    let open_area: f64 = mask
        .iter()
        .zip(surface.areas())
        .filter_map(|(&on, &a)| on.then_some(a))
        .sum();
    let structural: Vec<_> = surface
        .centroids()
        .iter()
        .zip(&mask)
        .map(|(x, &on)| if on { 1.0 + 0.5 * x[0] / 0.1 } else { 0.0 })
        .collect();
    let aperture: Vec<_> = mask
        .iter()
        .map(|&on| if on { 1.0 / open_area } else { 0.0 })
        .collect();
    let ids = [
        SourceId::StructuralMode("shell mode".into()),
        SourceId::ApertureVolumeFlow("opening".into()),
    ];
    let wall = LocalBoundary::try_new(
        medium,
        vec![
            LocalWall::Rigid,
            LocalWall::Mechanical(WallPin {
                surface_density: 0.1,
                stiffness_per_area: 5e4,
                resistance: 250.0,
            }),
        ],
        mask.iter().map(|&on| usize::from(!on)).collect(),
    )
    .unwrap();
    let model = observer(&surface, source, &ids, &[structural, aperture])
        .with_local_boundary(wall)
        .unwrap();
    let omega = 700.0;
    let response = model
        .prepare_response(omega, limits(&model), || false)
        .unwrap();
    assert_eq!(response.sources(), &ids);
    assert_eq!(response.source(), model.source());
    assert_eq!(response.harmonic().producer_hash(), model.operator_hash());
    assert!(response.acoustic_impedance().is_err());
    let rates = [C64::new(0.01, 0.007), C64::new(1e-4, -2e-4)];
    let drive: Vec<_> = (0..surface.areas().len())
        .map(|p| model.boundary.fields[0][p] * rates[0] + model.boundary.fields[1][p] * rates[1])
        .collect();
    let beta = model
        .local_boundary()
        .unwrap()
        .normalized_admittance(omega, &mut || false)
        .unwrap();
    let direct = model
        .exterior
        .solve_boundary_batch(
            &surface,
            omega,
            &[&drive],
            Admittance::Local {
                beta: &beta,
                declared_active: false,
            },
            Formulation::PlainCbie,
            &model.limits,
            || false,
        )
        .unwrap()
        .remove(0);
    let power = response.power_w(&rates, || false).unwrap();
    let scale = direct.driving_power_w.abs();
    close(power.driving_power_w, direct.driving_power_w, scale);
    close(power.absorbed_wall_power_w, direct.absorbed_power_w, scale);
    close(
        power.exterior_power_w,
        direct.radiation.boundary_work_w(),
        scale,
    );
    close(power.balance_residual_w, 0.0, scale);
    assert!(power.absorbed_wall_power_w > 0.0);
    let independent = response
        .power_w(&[rates[0], C64::ZERO], || false)
        .unwrap()
        .driving_power_w
        + response
            .power_w(&[C64::ZERO, rates[1]], || false)
            .unwrap()
            .driving_power_w;
    assert!((power.driving_power_w - independent).abs() > 1e-3 * scale);
    let mut efforts = [C64::ZERO; 2];
    response
        .resisting_efforts(&rates, &mut efforts, || false)
        .unwrap();
    close(
        0.5 * rates
            .iter()
            .zip(efforts)
            .map(|(&v, f)| (v.conj() * f).re)
            .sum::<f64>(),
        power.driving_power_w,
        scale,
    );
    let expected = direct
        .radiation
        .observe(&surface, &model.boundary.points, 0.01, || false)
        .unwrap();
    let mut actual = [HarmonicField::default()];
    response
        .harmonic()
        .evaluate_coherent(&rates, &mut actual, || false)
        .unwrap();
    assert!(
        (actual[0].pressure_pa - expected[0].pressure_pa).abs()
            < 2e-11 * expected[0].pressure_pa.abs()
    );
    for axis in 0..3 {
        let velocity =
            expected[0].gradient_pa_per_m[axis] * C64::new(0.0, -1.0 / (omega * medium.density));
        assert!(
            (actual[0].particle_velocity_m_s[axis] - velocity).abs()
                < 2e-11 * expected[0].pressure_pa.abs() / (medium.density * medium.sound_speed)
        );
    }
}

#[test]
fn preparation_limits_source_units_and_cancellation_refuse_without_new_solves() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let area: f64 = surface.areas().iter().sum();
    let model = observer(
        &surface,
        material(),
        &[SourceId::ApertureVolumeFlow("mouth".into())],
        &[vec![1.0 / area; surface.areas().len()]],
    );
    let admitted = limits(&model);
    let mut low_storage = admitted;
    low_storage.maximum_retained_entries -= 1;
    assert!(
        model
            .prepare_response(500.0, low_storage, || false)
            .is_err()
    );
    let mut low_work = admitted;
    low_work.maximum_panel_products -= 1;
    assert!(model.prepare_response(500.0, low_work, || false).is_err());
    for omega in [0.0, -1.0, f64::NAN, 3_001.0] {
        assert!(model.prepare_response(omega, admitted, || false).is_err());
    }
    assert!(matches!(
        model.prepare_response(500.0, admitted, || true),
        Err(MolecularObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
    let response = model.prepare_response(500.0, admitted, || false).unwrap();
    assert_eq!(response.acoustic_impedance().unwrap(), response.impedance());
    assert!(response.power_w(&[], || false).is_err());
    assert!(
        response
            .power_w(&[C64::new(f64::NAN, 0.0)], || false)
            .is_err()
    );
    assert!(matches!(
        response.power_w(&[C64::ONE], || true),
        Err(MolecularObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
    assert!(
        response
            .resisting_efforts(&[C64::ONE], &mut [], || false)
            .is_err()
    );
    let mut output = [C64::new(4.0, 5.0)];
    assert!(
        response
            .resisting_efforts(&[C64::ONE], &mut output, || true)
            .is_err()
    );
    assert_eq!(output, [C64::new(4.0, 5.0)]);
}
