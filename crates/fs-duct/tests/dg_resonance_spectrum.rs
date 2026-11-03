//! G0/G1/G3/G4: seed-free original DG modes, refinement and physical decay.
use fs_duct::dg::continuous::poles::{BeynOptions, DgPoleError, PoleOptions, ResonanceRegion};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, PortId};
use fs_duct::{Duct, Segment};
use fs_la::nonlinear_eig::beyn::BeynError;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_phs::WallPin;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn cylinder(
    elements: usize,
    order: usize,
    radius: f64,
    length: f64,
    inlet: DgBoundary,
    outlet: DgBoundary,
) -> DgDuct {
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        &air(),
        DgMesh { elements, order },
        inlet,
        outlet,
    )
    .unwrap()
}
fn options(rate: f64, probes: usize, nodes: usize) -> BeynOptions {
    BeynOptions {
        nodes,
        probe_columns: probes,
        rank_relative_tolerance: 1e-8,
        count_tolerance: 1e-5,
        quadrature_tolerance: 1e-5,
        contour_min_reciprocal_condition: 1e-10,
        linear_residual_tolerance: 1e-12,
        svd_tolerance: 1e-12,
        svd_max_sweeps: 128,
        eigen_max_sweeps: 256,
        max_dimension: 1024,
        max_dense_bytes: 64 * 1024 * 1024,
        max_work: 1_000_000_000_000,
        refinement: PoleOptions {
            spectral_scale: rate,
            step_tolerance: 1e-11,
            residual_tolerance: 1e-10,
            max_iterations: 24,
            max_dimension: 1024,
            max_dense_bytes: 32 * 1024 * 1024,
            max_work: 1_000_000_000,
        },
    }
}
fn pair_region(omega: f64) -> ResonanceRegion {
    ResonanceRegion {
        center: C64::new(0.0, -2.1 * omega),
        real_radius: omega,
        imaginary_radius: 1.65 * omega,
    }
}

#[test]
fn g1_two_acoustic_modes_are_discovered_without_seeds_and_converge_spatially() {
    let length = 0.5;
    let omega = core::f64::consts::PI * air().sound_speed / (2.0 * length);
    let region = pair_region(omega);
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let duct = cylinder(
            elements,
            2,
            0.01,
            length,
            DgBoundary::Port(PortId(7)),
            DgBoundary::IdealOpen,
        );
        let mut spectrum = duct
            .resonance_spectrum(region, &options(omega, 4, 128), || false)
            .unwrap();
        assert_eq!(spectrum.argument_count, 2);
        assert_eq!(spectrum.resolved_rank, 2);
        assert_eq!(spectrum.poles.len(), 2);
        assert_eq!(spectrum.quadrature_nodes, 128);
        spectrum.poles.sort_by(|a, b| {
            a.oscillation_frequency_rad_s
                .total_cmp(&b.oscillation_frequency_rad_s)
        });
        let mut rung = [0.0; 2];
        for (i, mode) in spectrum.poles.iter().enumerate() {
            let expected = (2 * i + 1) as f64 * omega;
            rung[i] = (mode.oscillation_frequency_rad_s / expected - 1.0).abs();
            assert!(mode.amplitude_decay_rate_per_s.abs() < 1e-10 * expected);
            assert!(mode.relative_residual < 1e-10);
            assert_eq!(mode.energy_coordinates.len(), duct.total_dofs());
            assert_eq!(mode.pressure_pa.len(), duct.field_dofs());
            assert_eq!(mode.volume_flow_m3_s.len(), duct.field_dofs());
            assert_eq!(mode.ports[0].id, PortId(7));
            assert!(mode.ports[0].pressure_pa.abs() > 0.0);
            assert_eq!(mode.wall_dissipation_w, 0.0);
            assert_eq!(mode.molecular_dissipation_w, 0.0);
            assert_eq!(mode.matched_absorption_w, 0.0);
            assert!((mode.stored_energy_j - 0.5).abs() < 2e-14);
        }
        errors.push(rung);
    }
    for mode in 0..2 {
        assert!(errors[1][mode] < 0.4 * errors[0][mode], "{errors:?}");
        assert!(errors[2][mode] < 0.4 * errors[1][mode], "{errors:?}");
        assert!(errors[2][mode] < 2e-3, "{errors:?}");
    }
}

#[test]
fn g1_contour_finds_material_relaxation_with_all_memory_and_original_heat() {
    let duct = cylinder(1, 1, 0.001, 0.05, DgBoundary::Closed, DgBoundary::Closed)
        .with_viscothermal(1, [10.0, 1000.0], 0.95)
        .unwrap();
    let wall = duct.viscothermal_model().unwrap();
    let term = wall.shunt_terms()[0];
    // Uniform sealed thermal mode: s = -(gain + rate), independently of DG.
    let expected = term.gain_per_s + term.rate_per_s;
    let region = ResonanceRegion {
        center: C64::from_re(-expected),
        real_radius: 0.04 * expected,
        imaginary_radius: 0.04 * expected,
    };
    let spectrum = duct
        .resonance_spectrum(region, &options(expected, 2, 64), || false)
        .unwrap();
    assert_eq!(spectrum.poles.len(), 1);
    assert_eq!(spectrum.argument_count, 1);
    let pole = &spectrum.poles[0];
    assert!((pole.laplace_pole_per_s + C64::from_re(expected)).abs() < 1e-9 * expected);
    assert_eq!(pole.oscillation_frequency_rad_s, 0.0);
    assert_eq!(pole.energy_coordinates.len(), duct.total_dofs());
    assert!(pole.energy_coordinates.len() > 2 * duct.field_dofs());
    assert!(pole.ports.is_empty());
    assert!(pole.auxiliary_stored_energy_j > 0.0);
    assert!(pole.wall_dissipation_w > 0.0);
    assert!(pole.decay_power_residual_w.abs() < 1e-10 * pole.wall_dissipation_w);
    let p = pole.pressure_pa[0];
    for value in &pole.pressure_pa {
        assert!((*value / p - C64::ONE).abs() < 1e-9);
    }
}

#[test]
fn g1_compliant_wall_contour_observes_independent_damped_oscillator_and_power() {
    let (radius, length) = (0.01, 0.3);
    let gas = air();
    let bulk = gas.density * gas.sound_speed.powi(2);
    let beta = 2.0 * bulk / radius;
    let wall = WallPin {
        surface_density: 0.2,
        stiffness_per_area: 1e7,
        resistance: 300.0,
    };
    let decay = wall.resistance / (2.0 * wall.surface_density);
    let omega = ((wall.stiffness_per_area + beta) / wall.surface_density - decay * decay).sqrt();
    let expected = C64::new(-decay, -omega);
    let duct = cylinder(1, 2, radius, length, DgBoundary::Closed, DgBoundary::Closed)
        .with_compliant_wall(wall)
        .unwrap();
    let region = ResonanceRegion {
        center: expected,
        real_radius: 0.01 * omega,
        imaginary_radius: 0.01 * omega,
    };
    let pole = duct
        .resonance_spectrum(region, &options(omega, 2, 64), || false)
        .unwrap()
        .poles
        .remove(0);
    assert!((pole.laplace_pole_per_s - expected).abs() < 1e-9 * omega);
    let displacement = pole.pressure_pa[0].scale(-1.0 / beta);
    let velocity = expected * displacement;
    let surface = core::f64::consts::TAU * radius * length;
    let heat = surface * wall.resistance * velocity.norm_sq();
    let stored = 0.5
        * surface
        * (wall.surface_density * velocity.norm_sq()
            + wall.stiffness_per_area * displacement.norm_sq());
    assert!((pole.wall_dissipation_w / heat - 1.0).abs() < 1e-9);
    assert!((pole.auxiliary_stored_energy_j / stored - 1.0).abs() < 1e-9);
    assert!(pole.decay_power_residual_w.abs() < 1e-9 * heat);
}

#[test]
fn g1_g3_prepared_spectrum_reuses_original_operator_with_live_quadrature() {
    let length = 0.5;
    let omega = core::f64::consts::PI * air().sound_speed / (2.0 * length);
    let duct = cylinder(
        1,
        2,
        0.01,
        length,
        DgBoundary::Port(PortId(7)),
        DgBoundary::IdealOpen,
    );
    let before = duct.continuous_operator().generator_per_s().to_vec();
    let prepared = duct
        .prepare_resonance_poles(omega, 1_000_000, || false)
        .unwrap();
    assert_eq!(prepared.dimension(), duct.total_dofs());
    let mut errors = Vec::new();
    for nodes in [16, 32, 64] {
        let mut controls = options(omega, 4, nodes);
        controls.count_tolerance = 0.2;
        controls.quadrature_tolerance = 0.2;
        controls.rank_relative_tolerance = 0.02;
        let spectrum = prepared
            .resonance_spectrum(pair_region(omega), &controls, || false)
            .unwrap();
        assert_eq!(spectrum.poles.len(), 2);
        let error = spectrum
            .unrefined_laplace_poles_per_s
            .iter()
            .map(|seed| {
                spectrum
                    .poles
                    .iter()
                    .map(|pole| (*seed - pole.laplace_pole_per_s).abs() / omega)
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0_f64, f64::max);
        errors.push(error);
    }
    assert!(
        errors[0] > 1e-12,
        "quadrature axis must be live: {errors:?}"
    );
    assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.2 * errors[1], "{errors:?}");
    let controls = options(omega, 4, 128);
    let one = prepared
        .resonance_spectrum(pair_region(omega), &controls, || false)
        .unwrap();
    let two = prepared
        .resonance_spectrum(pair_region(omega), &controls, || false)
        .unwrap();
    let direct = duct
        .resonance_spectrum(pair_region(omega), &controls, || false)
        .unwrap();
    assert_eq!(one, two);
    assert_eq!(one, direct);
    let first = &one.poles[0];
    let local = prepared
        .resonance_pole(
            first.laplace_pole_per_s,
            &first.energy_coordinates,
            &controls.refinement,
            || false,
        )
        .unwrap();
    assert!((local.laplace_pole_per_s - first.laplace_pole_per_s).abs() < 1e-10 * omega);
    assert_eq!(before, duct.continuous_operator().generator_per_s());
}

#[test]
fn g0_g4_missing_modes_resources_band_and_cancellation_refuse_transactionally() {
    let length = 0.5;
    let omega = core::f64::consts::PI * air().sound_speed / (2.0 * length);
    let duct = cylinder(
        1,
        2,
        0.01,
        length,
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
    );
    let before = duct.continuous_operator().generator_per_s().to_vec();
    let prepared = duct
        .prepare_resonance_poles(omega, 1_000_000, || false)
        .unwrap();
    assert!(matches!(
        prepared.resonance_spectrum(pair_region(omega), &options(omega, 1, 128), || false),
        Err(DgPoleError::Contour(BeynError::ProbeRankSaturated { .. }))
    ));
    let controls = options(omega, 4, 128);
    for mutated in [
        BeynOptions {
            max_dense_bytes: 1,
            ..controls
        },
        BeynOptions {
            max_work: 1,
            ..controls
        },
        BeynOptions {
            max_dimension: duct.total_dofs() - 1,
            ..controls
        },
        BeynOptions {
            probe_columns: duct.total_dofs() + 1,
            ..controls
        },
        BeynOptions {
            nodes: 0,
            ..controls
        },
    ] {
        assert!(
            prepared
                .resonance_spectrum(pair_region(omega), &mutated, || false)
                .is_err()
        );
    }
    assert!(duct.prepare_resonance_poles(omega, 1, || false).is_err());
    assert!(
        duct.prepare_resonance_poles(f64::NAN, 1_000_000, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(
        prepared
            .resonance_spectrum(pair_region(omega), &controls, || {
                polls += 1;
                polls > 30
            })
            .is_err()
    );
    assert!(polls > 30);
    assert_eq!(before, duct.continuous_operator().generator_per_s());
    // This duct's physical pole lies above its declared wall-law band. The
    // entire model is mathematically evaluable there, but must not be returned
    // as an admitted physical oscillatory mode.
    let narrow = cylinder(
        1,
        1,
        0.001,
        0.001,
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
    )
    .with_viscothermal(1, [10.0, 1000.0], 0.95)
    .unwrap();
    let fast = core::f64::consts::PI * air().sound_speed / (2.0 * 0.001);
    let outside = ResonanceRegion {
        center: C64::new(-1000.0, -fast),
        real_radius: 0.5 * fast,
        imaginary_radius: 0.5 * fast,
    };
    assert!(matches!(
        narrow.resonance_spectrum(outside, &options(fast, 3, 128), || false),
        Err(DgPoleError::Source(_))
    ));
}
