//! Source/material/receiver composition tests; spatial accuracy remains the
//! responsibility of the original BEM and material qualification ladders.
use super::*;
use crate::helmholtz::bc::{
    PorousRigidBoundary,
    incident::{self, IncidentField, Monopole},
};
use core::f64::consts::TAU;
use fs_material::fluid_acoustics::{FluidAcousticState, thermoviscous::ThermoviscousFluidState};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_phs::WallPin;

fn medium() -> Medium {
    Medium {
        density: 1.2,
        sound_speed: 340.0,
    }
}
fn source() -> PointSourceSpec {
    PointSourceSpec {
        id: "flow".into(),
        position_m: [0.02, 0.01, 0.09],
    }
}
fn receiver() -> ReceiverSpec {
    ReceiverSpec {
        id: "shadow-side".into(),
        position_m: [-0.03, 0.02, -0.08],
    }
}
fn scattering(surface: &SpherePanels) -> ScatteringObservation<'_> {
    ScatteringObservation::new(
        surface,
        medium(),
        &[source()],
        &[receiver()],
        TAU * 1600.0,
        1e-3,
        Options::default(),
        || false,
    )
    .unwrap()
}
fn pin() -> WallPin {
    WallPin {
        surface_density: 0.1,
        stiffness_per_area: 4e5,
        resistance: 200.0,
    }
}
fn porous_law() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap()
}
fn fluid_from_acoustic(acoustic: FluidAcousticState) -> ThermoviscousFluidState {
    ThermoviscousFluidState::from_heat_capacity_ratio(
        acoustic,
        FluidTransportState::try_new(acoustic.density_kg_m3(), 1005.0, 0.0257, 1.81e-5).unwrap(),
        1.4,
    )
    .unwrap()
}
fn fluid() -> ThermoviscousFluidState {
    fluid_from_acoustic(
        FluidAcousticState::try_new(medium().density, medium().sound_speed).unwrap(),
    )
}
fn porous() -> LocalWall {
    LocalWall::PorousRigid {
        law: porous_law(),
        fluid: fluid(),
        thickness_m: 0.03,
    }
}
fn uniform(surface: &SpherePanels, law: LocalWall) -> LocalBoundary {
    LocalBoundary::try_new(medium(), vec![law], vec![0; surface.areas().len()]).unwrap()
}
fn window() -> LowerBandWindow {
    LowerBandWindow::try_new(150.0, 400.0).unwrap()
}
fn declaration(owner: &LocalBoundaryObservation<'_>) -> TransferDeclaration {
    TransferDeclaration {
        producer_hash: owner.operator_hash(),
        sources: owner.sources().to_vec(),
        receivers: owner.receivers().to_vec(),
        fs_hz: 8000.0,
        geometric_delays_s: owner.geometric_delays_s(),
        pass_band_upper_hz: 800.0,
        taper_stop_hz: 1600.0,
        filter_latency_samples: 256,
        fft_lengths: vec![512, 1024, 2048, 4096],
        tail_energy_tolerance: 1e-10,
        dtft_relative_tolerance: 1e-4,
        producer_uncertainty: ProducerUncertainty::Unresolved {
            source_hash: owner.operator_hash(),
        },
        maximum_frequency_entries: 4 * MAX_OBSERVATION_FFT,
    }
}

#[test]
fn mechanical_and_porous_laws_keep_original_units_and_complex_phase() {
    let wall = pin();
    let boundary = LocalBoundary::try_new(
        medium(),
        vec![LocalWall::Rigid, LocalWall::Mechanical(wall), porous()],
        vec![2, 0, 1, 2],
    )
    .unwrap();
    let omega = TAU * 500.0;
    let beta = boundary.normalized_admittance(omega, || false).unwrap();
    // Independent mechanical impedance in exp(-i omega t); no conjugation,
    // area multiplication or guessed air impedance belongs in this conversion.
    let reactance = wall.stiffness_per_area / omega - omega * wall.surface_density;
    let denominator = wall.resistance.powi(2) + reactance.powi(2);
    let rho_c = medium().density * medium().sound_speed;
    let expected = C64::new(
        rho_c * wall.resistance / denominator,
        -rho_c * reactance / denominator,
    );
    assert!((beta[2] - expected).abs() < 1e-14);
    assert_eq!(beta[1], C64::ZERO);
    assert_eq!(beta[0], beta[3]);
    let law = porous_law();
    let gas = fluid();
    let original = PorousRigidBoundary::new(&law, &gas, omega, 0.03).unwrap();
    assert_eq!(beta[0], original.normalized_admittance());
    assert!(beta[0].re > 0.0 && beta[2].re > 0.0);
    assert_eq!(boundary.supported_band_hz(), (150.0, 1600.0));
    assert!(matches!(
        boundary.normalized_admittance(TAU * 100.0, || false),
        Err(TransferError::Material(_))
    ));
    assert!(matches!(
        boundary.normalized_admittance(omega, || true),
        Err(TransferError::Cancelled)
    ));
}

#[test]
fn complete_fluid_state_and_physical_patch_assignment_are_hash_bound() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let mut assignments = vec![0; surface.areas().len()];
    assignments[0] = 1;
    let first = LocalBoundary::try_new(
        medium(),
        vec![LocalWall::Rigid, LocalWall::Mechanical(pin())],
        assignments.clone(),
    )
    .unwrap();
    assignments.swap(0, 1);
    let second = LocalBoundary::try_new(
        medium(),
        vec![LocalWall::Rigid, LocalWall::Mechanical(pin())],
        assignments,
    )
    .unwrap();
    assert_ne!(first.content_hash(), second.content_hash());

    // An original K_s can differ while rounded sqrt(K_s/rho) is identical.
    // The porous compressibility consumes K_s, so rho/c alone cannot bind it.
    let first_bits = (medium().density * medium().sound_speed.powi(2)).to_bits();
    let pair = (0..32)
        .find_map(|offset| {
            let a = FluidAcousticState::from_adiabatic_bulk_modulus(
                medium().density,
                f64::from_bits(first_bits + offset),
            )
            .unwrap();
            let b = FluidAcousticState::from_adiabatic_bulk_modulus(
                medium().density,
                f64::from_bits(first_bits + offset + 1),
            )
            .unwrap();
            (a.sound_speed_m_s() == b.sound_speed_m_s()).then_some((a, b))
        })
        .expect("adjacent original bulk moduli sharing rounded sound speed");
    let shared = Medium {
        density: pair.0.density_kg_m3(),
        sound_speed: pair.0.sound_speed_m_s(),
    };
    let bind = |acoustic| {
        LocalBoundary::try_new(
            shared,
            vec![LocalWall::PorousRigid {
                law: porous_law(),
                fluid: fluid_from_acoustic(acoustic),
                thickness_m: 0.03,
            }],
            vec![0],
        )
        .unwrap()
    };
    assert_ne!(bind(pair.0).content_hash(), bind(pair.1).content_hash());
}

#[test]
fn point_source_scattering_keeps_rigid_parity_and_total_robin_field() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let omega = TAU * 700.0;
    let rigid = scattering(&surface)
        .frequency_sample(omega, || false)
        .unwrap();
    let unchanged = LocalBoundaryObservation::for_scattering(
        scattering(&surface),
        uniform(&surface, LocalWall::Rigid),
        None,
    )
    .unwrap()
    .physical_frequency_sample(omega, || false)
    .unwrap();
    assert_eq!(rigid, unchanged);
    let wall = uniform(&surface, LocalWall::Mechanical(pin()));
    let beta = wall.normalized_admittance(omega, || false).unwrap();
    let owner = LocalBoundaryObservation::for_scattering(scattering(&surface), wall, None).unwrap();
    let row = owner.physical_frequency_sample(omega, || false).unwrap();
    let monopoles = [Monopole {
        position_m: source().position_m,
        volume_velocity_m3_s: C64::ONE,
    }];
    let reference = incident::solve_scattering(
        &surface,
        omega / medium().sound_speed,
        medium(),
        Admittance::Local {
            beta: &beta,
            declared_active: false,
        },
        IncidentField::Monopoles(&monopoles),
        1e-3,
        helmholtz::Formulation::BurtonMiller,
    )
    .unwrap();
    assert!(reference.solution().absorbed_power_w > 0.0);
    assert!((row[0] - rigid[0]).abs() > 1e-3 * rigid[0].abs());
    let point = receiver().position_m;
    // Pressure parity uses the same vector quadrature owner. Scalar rows can
    // take another admitted refinement tree at their own error tolerance.
    let central = reference
        .total_pressure_velocity_at(&[point], 1e-3, Options::default())
        .unwrap();
    assert!((row[0] - central.scalar.pressure[0]).abs() < 1e-12 * row[0].abs());
    let h = 2e-5;
    let mut points = Vec::with_capacity(6);
    for axis in 0..3 {
        let mut plus = point;
        let mut minus = point;
        plus[axis] += h;
        minus[axis] -= h;
        points.extend([plus, minus]);
    }
    let pressure = reference
        .total_pressure_at(&points, 1e-3, Options::default())
        .unwrap();
    for axis in 0..3 {
        let gradient =
            (pressure.pressure[2 * axis] - pressure.pressure[1 + 2 * axis]).scale(0.5 / h);
        assert!((row[axis + 1] - gradient).abs() < 2e-5 * gradient.abs());
    }
}

#[test]
fn finite_apertures_keep_prescribed_flow_while_other_patches_absorb() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let mut traces = vec![vec![0.0; surface.areas().len()]; 2];
    for index in 0..2 {
        traces[index][index] = 1.0 / surface.areas()[index];
    }
    let sources = [
        SourceId::ApertureVolumeFlow("a".into()),
        SourceId::ApertureVolumeFlow("b".into()),
    ];
    let mut patches = vec![1; surface.areas().len()];
    patches[..2].fill(0);
    let wall = LocalBoundary::try_new(
        medium(),
        vec![LocalWall::Rigid, LocalWall::Mechanical(pin())],
        patches,
    )
    .unwrap();
    let original =
        BemObservation::new(&surface, medium(), &sources, &traces, &[receiver()], 1e-3).unwrap();
    let owner = LocalBoundaryObservation::for_surface(
        original,
        wall.clone(),
        TAU * 1600.0,
        Options::default(),
        None,
        || false,
    )
    .unwrap();
    let omega = TAU * 500.0;
    let observed = owner.physical_frequency_sample(omega, || false).unwrap();
    assert_eq!(observed.len(), 8);
    let beta = wall.normalized_admittance(omega, || false).unwrap();
    let prepared = PreparedBoundary::new(
        &surface,
        omega / medium().sound_speed,
        medium(),
        Admittance::Local {
            beta: &beta,
            declared_active: false,
        },
        helmholtz::Formulation::BurtonMiller,
    )
    .unwrap();
    let points = [receiver().position_m];
    let geometry = Geometry::new(&surface, &points, 1e-3).unwrap();
    let rows = geometry
        .prepare_velocity(omega / medium().sound_speed, medium(), Options::default())
        .unwrap();
    for (index, trace) in traces.iter().enumerate() {
        let field: Vec<_> = trace.iter().map(|&value| C64::from_re(value)).collect();
        let response = prepared.solve(&field).unwrap();
        assert_eq!(response.radiation.velocity[index], field[index]);
        assert!(
            (response.radiation.velocity[index].scale(surface.areas()[index]) - C64::ONE).abs()
                < 1e-14
        );
        assert!(response.absorbed_power_w > 0.0);
        assert!(
            response.radiation.velocity[2..]
                .iter()
                .any(|&value| value != C64::ZERO)
        );
        let received = rows.evaluate_velocity(&response.radiation).unwrap();
        assert_eq!(observed[4 * index], received.scalar.pressure[0]);
        for axis in 0..3 {
            let gradient =
                C64::new(0.0, omega * medium().density) * received.particle_velocity_m_s[0][axis];
            assert_eq!(observed[4 * index + axis + 1], gradient);
        }
    }
    let original =
        BemObservation::new(&surface, medium(), &sources, &traces, &[receiver()], 1e-3).unwrap();
    assert!(matches!(
        LocalBoundaryObservation::for_surface(
            original,
            uniform(&surface, LocalWall::Mechanical(pin())),
            TAU * 1600.0,
            Options::default(),
            None,
            || false,
        ),
        Err(TransferError::Invalid(_))
    ));
    let ruled = BemObservation::new_ruled_burton_miller(
        &surface,
        medium(),
        &sources,
        &traces,
        &[receiver()],
        1e-3,
        helmholtz::PanelRule::Triangle { max_pairs: 4096 },
    )
    .unwrap();
    assert!(matches!(
        LocalBoundaryObservation::for_surface(
            ruled,
            wall,
            TAU * 1600.0,
            Options::default(),
            None,
            || false,
        ),
        Err(TransferError::Invalid(_))
    ));
}

#[test]
fn unsupported_material_bands_and_lossless_poles_refuse_before_solving() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    assert!(matches!(
        LocalBoundaryObservation::for_scattering(
            scattering(&surface),
            uniform(&surface, porous()),
            None,
        ),
        Err(TransferError::Invalid(_))
    ));
    let owner = LocalBoundaryObservation::for_scattering(
        scattering(&surface),
        uniform(&surface, porous()),
        Some(window()),
    )
    .unwrap();
    assert!(matches!(
        owner.physical_frequency_sample(TAU * 100.0, || false),
        Err(TransferError::Material(_))
    ));
    assert_eq!(
        owner
            .band_limited_frequency_sample(TAU * 100.0, || false)
            .unwrap(),
        vec![C64::ZERO; 4]
    );
    assert!(matches!(
        owner.band_limited_frequency_sample(0.0, || true),
        Err(TransferError::Cancelled)
    ));
    assert!(
        owner
            .band_limited_frequency_sample(TAU * 1601.0, || false)
            .is_err()
    );
    let lossless = uniform(
        &surface,
        LocalWall::Mechanical(WallPin {
            surface_density: 1.0,
            stiffness_per_area: 1e6,
            resistance: 0.0,
        }),
    );
    assert!(matches!(
        lossless.normalized_admittance(1000.0, || false),
        Err(TransferError::Invalid(_))
    ));
    assert!(matches!(
        LocalBoundaryObservation::for_scattering(scattering(&surface), lossless.clone(), None,),
        Err(TransferError::Invalid(_))
    ));
    assert!(
        LocalBoundaryObservation::for_scattering(
            scattering(&surface),
            lossless,
            Some(LowerBandWindow::try_new(200.0, 400.0).unwrap()),
        )
        .is_ok()
    );
    let wrong_fluid = fluid_from_acoustic(FluidAcousticState::try_new(1.1, 340.0).unwrap());
    assert!(matches!(
        LocalBoundary::try_new(
            medium(),
            vec![LocalWall::PorousRigid {
                law: porous_law(),
                fluid: wrong_fluid,
                thickness_m: 0.03,
            }],
            vec![0]
        ),
        Err(TransferError::Invalid(_))
    ));
}

#[test]
fn lower_window_requires_live_passband_bins_and_binds_artifact_identity() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let owner = LocalBoundaryObservation::for_scattering(
        scattering(&surface),
        uniform(&surface, porous()),
        Some(window()),
    )
    .unwrap();
    let other = LocalBoundaryObservation::for_scattering(
        scattering(&surface),
        uniform(&surface, porous()),
        Some(LowerBandWindow::try_new(150.0, 450.0).unwrap()),
    )
    .unwrap();
    assert_ne!(owner.operator_hash(), other.operator_hash());
    let window = window();
    assert_eq!(window.weight(150.0), 0.0);
    assert_eq!(window.weight(400.0), 1.0);
    assert_eq!(window.weight(275.0), 0.5);
    for frequency in [151.0, 170.0, 275.0, 390.0, 399.0] {
        assert!(
            (window.weight(frequency) + observation_taper(frequency, 150.0, 400.0) - 1.0).abs()
                < 1e-14
        );
    }
    let mut narrow = declaration(&owner);
    narrow.pass_band_upper_hz = 410.0;
    assert!(narrow.admit_work().is_ok());
    assert!(matches!(
        owner.bake(&narrow, || false),
        Err(TransferError::Invalid(_))
    ));
    let d = declaration(&owner);
    assert!(matches!(
        other.bake(&d, || false),
        Err(TransferError::Invalid(_))
    ));
    assert!(matches!(
        owner.bake(&d, || true),
        Err(TransferError::Cancelled)
    ));
}

#[test]
fn porous_diffraction_bakes_to_explicit_band_limited_pressure_playback() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let wall = uniform(&surface, porous());
    let hash = wall.content_hash();
    let owner =
        LocalBoundaryObservation::for_scattering(scattering(&surface), wall, Some(window()))
            .unwrap();
    let d = declaration(&owner);
    let bake = owner.bake(&d, || false).unwrap();
    assert_eq!(bake.lower_band_window(), Some(window()));
    assert_eq!(bake.boundary_hash(), hash);
    let transfer = bake.transfer();
    assert_eq!(transfer.declaration_hash, d.declaration_hash());
    assert_eq!(transfer.producer_uncertainty, d.producer_uncertainty);
    assert_eq!(transfer.ladder.len(), 4);
    let omega = TAU * 500.0;
    let sample = owner.physical_frequency_sample(omega, || false).unwrap();
    let latency = phase(-omega * d.filter_latency_samples as f64 / d.fs_hz);
    for (kernel, sample) in transfer.kernels.iter().zip(sample) {
        let reference = sample.conj() * latency;
        let measured = kernel.transfer_at(omega / d.fs_hz).unwrap();
        assert!((measured - reference).abs() < 2e-4 * reference.abs().max(1.0));
    }
    let prefix = (transfer.geometric_delays_s[0] * transfer.fs_hz).floor() as usize;
    let mut observer = bake.observer().unwrap();
    let mut peak = 0.0_f64;
    for time in 0..512 {
        let output = observer
            .step(&[if time == 0 { 1e-6 } else { 0.0 }])
            .unwrap();
        assert_eq!(output.len(), 4);
        assert!(output.iter().all(|value| value.is_finite()));
        if time < prefix {
            assert_eq!(output, vec![0.0; 4]);
        }
        peak = peak.max(output[0].abs());
    }
    assert!(peak > 1e-4);
}
