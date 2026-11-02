use super::*;
use core::f64::consts::{PI, TAU};
use fs_material::gas::GasSpec;
use fs_math::{c64::C64, det};
use fs_vfit::mimo::{MimoModel, load::MimoLoadLimits};
use fs_vfit::model::PoleTerm;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn budget() -> RationalModalLineBudget {
    RationalModalLineBudget {
        linear: LinearSplitBudget {
            max_modes: 124,
            max_omega_dt: 1.0,
        },
        exact: DenseExponentialBudget {
            max_dimension: 256,
            max_dense_bytes: 64 << 20,
            max_multiply_terms: 20_000_000_000,
        },
        relative_tolerance: 1e-8,
    }
}
fn reference() -> f64 {
    gas().density * gas().sound_speed / (PI * 0.025_f64.powi(2))
}
fn model(ports: usize, direct: bool, pair: bool) -> MimoModel {
    let z0 = reference();
    let mutual = if ports == 1 {
        vec![0.4]
    } else {
        vec![0.8, 0.28, 0.28, 0.6]
    };
    let mut poles = vec![PoleTerm::Real {
        pole: -3500.0,
        residue: 0.0,
    }];
    let mut residues = vec![
        mutual
            .iter()
            .map(|r| C64::from_re(r * z0 * 3500.0))
            .collect(),
    ];
    if pair {
        // Positive-real damped resonator: 2R*s / ((s+alpha)^2+beta^2).
        poles.push(PoleTerm::Pair {
            pole: C64::new(-400.0, 12000.0),
            residue: C64::ZERO,
        });
        residues.push(
            mutual
                .iter()
                .map(|r| C64::new(r * z0 * 400.0, r * z0 * 400.0 * 400.0 / 12000.0))
                .collect(),
        );
    }
    MimoModel {
        ports,
        poles,
        residues,
        d: mutual
            .iter()
            .map(|v| if direct { v * z0 * 0.2 } else { 0.0 })
            .collect(),
        e: vec![0.0; ports * ports],
        reciprocal: true,
    }
}
fn load(model: &MimoModel, dt: f64) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        model,
        dt,
        MimoLoadLimits {
            band_rad_s: [1.0, 4.0e6],
            sampled_passivity_points: 33,
            max_states: 32,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}
fn line(model: &MimoModel, dt: f64, length: f64, elements: usize) -> RationalModalLine {
    let profile = RadiusProfile::exponential(0.025, 0.0, length).unwrap();
    RationalModalLine::new(
        &profile,
        &gas(),
        elements,
        ModalRationalLoad::new(load(model, dt), 0.025, gas()).unwrap(),
        budget(),
        &CancelGate::new(),
    )
    .unwrap()
}
fn close(a: f64, b: f64, scale: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * scale.max(1e-28),
        "{a:e} != {b:e}, scale={scale:e}"
    );
}
fn snapshot(line: &RationalModalLine) -> Vec<f64> {
    let mut values = line.state().to_vec();
    values.extend(line.incoming());
    for port in line.port_accounts() {
        values.extend([
            port.boundary_work,
            port.incident_energy,
            port.outgoing_energy,
            port.outgoing_wave,
        ]);
    }
    for port in line.terminal_accounts() {
        values.extend([
            port.mean_duct_pressure_pa,
            port.mean_load_pressure_pa,
            port.mean_outward_flow_m3_s,
            port.duct_work_j,
            port.load_work_j,
        ]);
    }
    for state in line.load().pole_states() {
        values.extend([state.re, state.im]);
    }
    values
}

#[test]
fn complete_zero_direct_mutual_mouth_evolves_original_history_and_independent_work() {
    let source = model(2, false, true);
    let mut region = line(&source, 1.0 / 96_000.0, 0.035, 3);
    let gate = CancelGate::new();
    assert_eq!(region.state().len(), 16 + source.order() * 2);
    assert_eq!(region.load().source_model(), &source);
    let (continuous_class, sampled_class) = (
        region.load().reports().0.class,
        region.load().reports().1.class,
    );
    let (mut incoming, mut outgoing, mut terminal) = (0.0, 0.0, 0.0);
    let (mut mutual_live, mut returned_work, mut varying_wave) = (false, false, false);
    let mut stored_at_release = 0.0;
    for step in 0..480 {
        // Only the higher, zero-net-flow mode is driven. A cylindrical bore has
        // no axial mode conversion, so any plane reaction comes from the mouth.
        let wave = if step < 48 {
            [0.0, 70.0 * det::sin(0.14 * (step as f64 + 0.5))]
        } else {
            [0.0; 2]
        };
        let mut free = [0.0; 2];
        let mut direct = [0.0; 4];
        region
            .reflection_relation_into(&mut free, &mut direct)
            .unwrap();
        let expected = [
            free[0] + direct[0] * wave[0] + direct[1] * wave[1],
            free[1] + direct[2] * wave[0] + direct[3] * wave[1],
        ];
        let record = region.push(&wave, &gate).unwrap();
        for k in 0..2 {
            close(region.incoming()[k], expected[k], 100.0, 2e-9);
            let port = region.port_accounts()[k];
            close(
                port.incident_energy,
                region.duration() * wave[k].powi(2) / reference(),
                port.incident_energy.abs(),
                2e-11,
            );
            varying_wave |= port.outgoing_energy
                - region.duration() * port.outgoing_wave.powi(2) / reference()
                > 1e-13;
            let mouth = region.terminal_accounts()[k];
            close(
                mouth.mean_duct_pressure_pa,
                mouth.mean_load_pressure_pa,
                100.0,
                2e-9,
            );
            close(
                mouth.duct_work_j,
                mouth.load_work_j,
                record.incident_energy_j
                    + record.reflected_energy_j
                    + region.stored_duct_energy_j(),
                2e-8,
            );
            returned_work |= mouth.load_work_j < -1e-13;
        }
        mutual_live |= region.incoming()[0].abs() > 1e-3;
        incoming += record.incident_energy_j;
        outgoing += record.reflected_energy_j;
        terminal += record.terminal_work_j;
        close(
            region.stored_duct_energy_j(),
            incoming - outgoing - terminal,
            incoming + outgoing + terminal.abs(),
            2e-8,
        );
        assert_eq!(
            region.load_state(),
            region.load().continuous_operator().real_coordinates()
        );
        if step == 48 {
            stored_at_release = region.stored_duct_energy_j();
        }
    }
    assert!(mutual_live && returned_work && varying_wave);
    assert!(region.stored_duct_energy_j() < 0.1 * stored_at_release);
    assert!(region.load().pole_states().iter().any(|z| z.im != 0.0));
    assert_eq!(region.load().source_model(), &source);
    assert_eq!(region.load().reports().0.class, continuous_class);
    assert_eq!(region.load().reports().1.class, sampled_class);
    assert_eq!(region.duration(), 1.0 / 96_000.0);
}

#[test]
fn held_wave_subdivision_preserves_every_original_coordinate_and_signed_exchange() {
    let source = model(2, true, true);
    let dt = 1.0 / 24_000.0;
    let mut coarse = line(&source, dt, 0.035, 3);
    let mut half = line(&source, dt / 2.0, 0.035, 3);
    let mut quarter = line(&source, dt / 4.0, 0.035, 3);
    let gate = CancelGate::new();
    for step in 0..96 {
        let wave = if step < 24 {
            [
                50.0 * det::sin(0.21 * step as f64),
                30.0 * det::cos(0.17 * step as f64),
            ]
        } else {
            [0.0; 2]
        };
        let accepted = coarse.push(&wave, &gate).unwrap();
        for (refined, divisions) in [(&mut half, 2), (&mut quarter, 4)] {
            let mut work = 0.0;
            let mut reflected = 0.0;
            let mut means = [0.0; 2];
            for _ in 0..divisions {
                let record = refined.push(&wave, &gate).unwrap();
                work += record.terminal_work_j;
                reflected += record.reflected_energy_j;
                for k in 0..2 {
                    means[k] += refined.incoming()[k] / divisions as f64;
                }
            }
            let q = coarse.source().coordinates();
            for range in [0..q, q..2 * q, 2 * q..coarse.state().len()] {
                let scale = coarse.state()[range.clone()]
                    .iter()
                    .chain(&refined.state()[range.clone()])
                    .map(|v| v.abs())
                    .fold(0.0, f64::max);
                for i in range {
                    close(coarse.state()[i], refined.state()[i], scale, 3e-8);
                }
            }
            let scale = coarse.stored_duct_energy_j()
                + accepted.incident_energy_j
                + accepted.reflected_energy_j;
            close(work, accepted.terminal_work_j, scale, 3e-8);
            close(reflected, accepted.reflected_energy_j, scale, 3e-8);
            for k in 0..2 {
                close(means[k], coarse.incoming()[k], 100.0, 3e-9);
            }
        }
    }
}

// Independent weak harmonic equation, without the new joint generator or
// exact stepper: (K-w²M-iw B_in B_in^T/Z0-iw B_out B_out^T/ZL) q=2 B_in/Z0.
fn weak_reflection(line: &RationalModalLine, omega: f64) -> C64 {
    let source = line.source();
    let n = source.coordinates();
    let s = C64::new(0.0, -omega);
    let z = line.load().source_model().eval(s)[0];
    let mut pencil: Vec<_> = source
        .stiffness()
        .iter()
        .zip(source.mass())
        .map(|(k, m)| C64::from_re(k - omega * omega * m))
        .collect();
    pencil[0] = pencil[0] + s.scale(1.0 / source.inlet_reference_impedance());
    pencil[n * n - 1] = pencil[n * n - 1] + s / z;
    let mut response = vec![C64::ZERO; n];
    response[0] = C64::from_re(2.0 / source.inlet_reference_impedance());
    fs_la::eigen_complex::lu_complex(&pencil, n)
        .unwrap()
        .solve(&mut response);
    s * response[0] - C64::ONE
}

#[test]
fn reflected_harmonic_converges_to_original_weak_equation_and_analytic_cylinder() {
    let source = model(1, false, false);
    let omega = TAU * 1100.0;
    let length = 0.12;
    let terminal = source.eval(C64::new(0.0, -omega))[0];
    let phase = 2.0 * omega * length / gas().sound_speed;
    let exact = (terminal - C64::from_re(reference())) / (terminal + C64::from_re(reference()))
        * C64::new(det::cos(phase), det::sin(phase));
    let gate = CancelGate::new();
    let (mut spatial, mut temporal) = (Vec::new(), Vec::new());
    for (elements, samples) in [(3, 64), (6, 128), (12, 256)] {
        let dt = TAU / (omega * samples as f64);
        let mut region = line(&source, dt, length, elements);
        let weak = weak_reflection(&region, omega);
        let mut measured = C64::ZERO;
        for step in 0..12 * samples {
            let phase = omega * dt * (step as f64 + 0.5);
            region.push(&[det::cos(phase)], &gate).unwrap();
            if step >= 8 * samples {
                measured = measured
                    + C64::new(det::cos(phase), det::sin(phase))
                        .scale(2.0 * region.incoming()[0] / (4 * samples) as f64);
            }
        }
        spatial.push((weak - exact).abs());
        temporal.push((measured - weak).abs());
    }
    assert!(
        spatial[1] < 0.32 * spatial[0] && spatial[2] < 0.32 * spatial[1],
        "{spatial:?}"
    );
    assert!(
        temporal[1] < 0.4 * temporal[0] && temporal[2] < 0.4 * temporal[1],
        "{temporal:?}"
    );
    assert!(
        spatial[2] < 0.012 && temporal[2] < 0.0002,
        "spatial={spatial:?}; temporal={temporal:?}"
    );
    eprintln!(
        "rational modal cylinder spatial errors={spatial:?}; held-wave harmonic errors={temporal:?}"
    );
}

#[test]
fn original_history_refusals_and_failed_pushes_are_transactional() {
    let source = model(2, true, true);
    let dt = 1.0 / 96_000.0;
    let profile = RadiusProfile::exponential(0.025, 0.0, 0.035).unwrap();
    let gate = CancelGate::new();
    let mut original = load(&source, dt);
    original
        .step_mean_into(&[1e-5, -2e-5], &mut [0.0; 2])
        .unwrap();
    let accepted = original.continuous_operator().real_coordinates();
    let mut region = RationalModalLine::new(
        &profile,
        &gas(),
        3,
        ModalRationalLoad::new(original, 0.025, gas()).unwrap(),
        budget(),
        &gate,
    )
    .unwrap();
    assert_eq!(region.load_state(), accepted);
    region.push(&[10.0, -30.0], &gate).unwrap();
    let before = snapshot(&region);
    for bad in [&[f64::MAX, f64::MAX][..], &[f64::NAN, 0.0][..], &[1.0][..]] {
        assert!(region.push(bad, &gate).is_err());
        assert_eq!(snapshot(&region), before);
    }
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        region.push(&[0.0; 2], &cancelled),
        Err(RationalModalError::Cancelled)
    ));
    assert_eq!(snapshot(&region), before);
    let mut free = [f64::MAX];
    let mut direct = [f64::MAX; 4];
    assert!(
        region
            .reflection_relation_into(&mut free, &mut direct)
            .is_err()
    );
    assert_eq!(free, [f64::MAX]);
    assert_eq!(direct, [f64::MAX; 4]);

    let mut zero = load(&model(2, false, false), dt);
    zero.step_mean_into(&[1e-5, 0.0], &mut [0.0; 2]).unwrap();
    assert!(
        RationalModalLine::new(
            &profile,
            &gas(),
            3,
            ModalRationalLoad::new(zero, 0.025, gas()).unwrap(),
            budget(),
            &gate
        )
        .is_err()
    );
    let mut mixed = model(2, false, false);
    mixed.d[0] = reference();
    let error = RationalModalLine::new(
        &profile,
        &gas(),
        3,
        ModalRationalLoad::new(load(&mixed, dt), 0.025, gas()).unwrap(),
        budget(),
        &gate,
    )
    .err()
    .unwrap();
    assert!(error.to_string().contains("mixed-rank"), "{error}");
    assert!(
        RationalModalLine::new(
            &profile,
            &gas(),
            3,
            ModalRationalLoad::new(load(&source, dt), 0.026, gas()).unwrap(),
            budget(),
            &gate
        )
        .is_err()
    );
    let mut tiny = budget();
    tiny.exact.max_multiply_terms = 1;
    assert!(
        RationalModalLine::new(
            &profile,
            &gas(),
            3,
            ModalRationalLoad::new(load(&source, dt), 0.025, gas()).unwrap(),
            tiny,
            &gate
        )
        .is_err()
    );
    // This fits the initial usize dimension arithmetic but would overflow an
    // unchecked cubic work estimate. Refuse before dense-byte/work products.
    assert!(
        RationalModalLine::new(
            &profile,
            &gas(),
            usize::MAX / 8,
            ModalRationalLoad::new(load(&source, dt), 0.025, gas()).unwrap(),
            budget(),
            &gate,
        )
        .is_err()
    );
}

#[test]
fn original_geometric_bem_realization_drives_the_complete_continuous_modal_duct() {
    use crate::bernoulli_aperture::radiation::environment::multiport::{
        MimoRadiationLimits, bake::BaffledMimoBakeOptions,
    };
    use crate::modal_mouth::{
        BemModalMouth, CircularMouthPatch, ModalMouthBemBudget,
        realization::ModalMouthRealizationOptions,
    };
    use fs_bem::{helmholtz::aperture::ApertureRoute, panel3d::SpherePanels};
    use fs_vfit::{causality::CausalityOptions, passivity::CertificateClass};

    // Coarse physical closed cylinder with two cap annuli. This qualifies the
    // native source-to-runtime connection, not continuum or out-of-band accuracy.
    let radius = 0.02;
    let bell =
        fs_bem::radiation_bake::lathe_profile(&[(0.0, radius), (0.01, radius), (0.02, radius)], 8)
            .unwrap();
    let patch = CircularMouthPatch {
        centre_m: [0.02, 0.0, 0.0],
        outward_normal: [1.0, 0.0, 0.0],
        radius_m: radius,
        panels: (bell.mouth_start..bell.triangles.len()).collect(),
    };
    let surface = SpherePanels::from_triangles(bell.triangles).unwrap();
    let gate = CancelGate::new();
    let source = BemModalMouth::new(
        &surface,
        patch,
        gas(),
        ModeCount(2),
        ModalMouthBemBudget {
            maximum_panels: 200,
            maximum_trace_entries: 400,
            maximum_matrix_entries: 2000,
            maximum_projection_terms: 1_000_000,
            maximum_basis_gram_defect: 0.45,
        },
        &gate,
    )
    .unwrap();
    let realized = source
        .realize(
            &ApertureRoute::Galerkin {
                budget: fs_bem::reciprocal_radiation::Budget {
                    maximum_unknowns: 400,
                    maximum_cache_bytes: 16 << 20,
                    maximum_work: 4_000_000_000,
                    angular_refinement: 0,
                },
            },
            ModalMouthRealizationOptions {
                mechanics_rate_hz: 6000,
                source_band_hz: [1.0, 3000.0],
                operator_band_hz: [30.0, 3000.0],
                fit: BaffledMimoBakeOptions {
                    frequency_samples: 32,
                    orders: [2, 4, 6],
                    fit_iterations: 12,
                    maximum_relative_error: 0.03,
                    max_fit_memory_bytes: 64 << 20,
                    causality: CausalityOptions {
                        subtractions: 2,
                        tail_bound_ratio: 100.0,
                        relative_data_tolerance: 1e-8,
                    },
                },
                load: MimoRadiationLimits {
                    max_states: 32,
                    max_memory_bytes: 1 << 20,
                    sampled_passivity_points: 257,
                },
                maximum_reciprocity_defect: 1e-10,
            },
            &gate,
        )
        .unwrap();
    let original = realized.load().source_model().clone();
    let profile = RadiusProfile::exponential(radius, 0.0, 0.06).unwrap();
    let mut line = RationalModalLine::new(
        &profile,
        &gas(),
        3,
        ModalRationalLoad::from_bem(realized),
        budget(),
        &gate,
    )
    .unwrap();
    assert_eq!(line.load().source_model(), &original);
    assert_eq!(line.state().len(), 16 + 2 * original.order());
    assert_eq!(line.duration(), 1.0 / 6000.0);
    let record = line.boundary().radiation_record().unwrap();
    assert_eq!(record.training().loads().len(), 32);
    assert_eq!(record.held_out().loads().len(), 31);
    assert_eq!(record.causality().loads().len(), 32);
    assert_eq!(record.options().source_band_hz, [1.0, 3000.0]);
    assert_eq!(record.report().sampled.class, CertificateClass::GridOnly);
    let mut total_in = 0.0;
    let mut total_out = 0.0;
    let mut total_terminal = 0.0;
    let mut mutual_live = false;
    for step in 0..80 {
        let incident = if step < 24 {
            [
                0.0,
                20.0 * det::sin(TAU * 700.0 * line.duration() * (step as f64 + 0.5)),
            ]
        } else {
            [0.0; 2]
        };
        let exchange = line.push(&incident, &gate).unwrap();
        total_in += exchange.incident_energy_j;
        total_out += exchange.reflected_energy_j;
        total_terminal += exchange.terminal_work_j;
        mutual_live |= line.incoming()[0].abs() > 1e-4;
        close(
            line.stored_duct_energy_j(),
            total_in - total_out - total_terminal,
            total_in + total_out + total_terminal.abs(),
            3e-8,
        );
    }
    assert!(mutual_live);
    assert_eq!(line.load().source_model(), &original);
    assert_eq!(
        line.load_state(),
        line.load().continuous_operator().real_coordinates()
    );
    assert!(line.load().pole_states().iter().any(|z| z.abs() > 0.0));
    assert_eq!(
        line.boundary().radiation_record().unwrap().modes(),
        ModeCount(2)
    );
}
