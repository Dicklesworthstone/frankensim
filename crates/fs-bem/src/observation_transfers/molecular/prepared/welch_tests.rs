//! Retained molecular input/exterior/wall work under measured source coherence.
use super::*;
use crate::helmholtz::Medium;
use crate::observation_transfers::ReceiverSpec;
use fs_blake3::ContentHash;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};

const N: usize = 32;
const BIN: usize = 3;
const DT: f64 = 1.0 / 4096.0;
const WORK: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 66,
};

fn material() -> RelaxationState {
    // Existing declared caloric-reservoir fixture, not measured gas parameters.
    let text = "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500";
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}

fn spectrum(segments: &[[C64; 2]], window: WelchWindow) -> WelchCrossSpectrum {
    let rows: Vec<_> = segments
        .iter()
        .flat_map(|values| {
            (0..N).map(move |j| {
                let theta = core::f64::consts::TAU * BIN as f64 * j as f64 / N as f64;
                values.map(|a| a.re * theta.cos() + a.im * theta.sin())
            })
        })
        .collect();
    WelchCrossSpectrum::estimate_interleaved(
        &rows,
        DT,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 2,
            maximum_samples_per_channel: 3 * N,
            maximum_transforms: 6,
            maximum_workspace_bytes: 100_000,
        },
        || false,
    )
    .unwrap()
}

fn segments() -> [[C64; 2]; 3] {
    [
        [C64::new(0.003, -0.002), C64::new(2e-5, 1e-5)],
        [C64::new(-0.001, 0.004), C64::new(-1e-5, 2e-5)],
        [C64::new(0.002, 0.001), C64::new(3e-5, -1e-5)],
    ]
}

fn owner(spectrum: &WelchCrossSpectrum) -> PreparedMolecularResponse {
    let source = material();
    let medium = Medium {
        density: source.gas().density,
        sound_speed: source.gas().sound_speed,
    };
    let ids = [
        SourceId::StructuralMode("surface".into()),
        SourceId::ApertureVolumeFlow("opening".into()),
    ];
    let harmonic = PreparedHarmonicObservation::from_samples(
        core::f64::consts::TAU * spectrum.frequency_hz(BIN).unwrap(),
        medium,
        &ids,
        &[ReceiverSpec {
            id: "unused receiver".into(),
            position_m: [0.0, 0.0, 1.0],
        }],
        ContentHash([43; 32]),
        vec![C64::ZERO; 8],
    )
    .unwrap();
    // Raw retained matrices for algebraic query checks; no BEM solve or
    // passivity/continuum accuracy is inferred from this fixture.
    let impedance = vec![
        C64::new(800.0, 300.0),
        C64::new(10_000.0, 2_000.0),
        C64::new(8_000.0, -1_000.0),
        C64::new(4e6, 2e6),
    ];
    PreparedMolecularResponse {
        source,
        harmonic,
        boundary_solutions: vec![],
        receiver_quadrature: None,
        exterior_work: impedance.iter().map(|z| z.scale(0.7)).collect(),
        wall_work: impedance.iter().map(|z| z.scale(0.3)).collect(),
        impedance,
    }
}

fn oracle(matrix: &[C64], rates: [C64; 2]) -> f64 {
    let mut power = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            power += 0.5 * (rates[i].conj() * matrix[2 * i + j] * rates[j]).re;
        }
    }
    power
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-20),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn measured_molecular_work_keeps_phase_mixed_units_and_separate_power_channels() {
    let inputs = segments();
    for (window, fraction) in [
        (WelchWindow::Rectangular, 1.0),
        (WelchWindow::PeriodicHann, 2.0 / 3.0),
    ] {
        let spectrum = spectrum(&inputs, window);
        let owner = owner(&spectrum);
        let expected = [
            owner.impedance(),
            owner.exterior_work_matrix(),
            owner.wall_work_matrix(),
        ]
        .map(|matrix| inputs.iter().map(|&q| oracle(matrix, q)).sum::<f64>() * fraction / 3.0);
        let mut scratch = [C64::ZERO; 2];
        let actual = owner
            .power_welch_bin(&spectrum, BIN, owner.sources(), WORK, &mut scratch, || {
                false
            })
            .unwrap();
        close(actual.driving_power_w, expected[0]);
        close(actual.exterior_power_w, expected[1]);
        close(actual.absorbed_wall_power_w, expected[2]);
        assert!(actual.balance_residual_w.abs() <= 3e-12 * actual.driving_power_w.abs());
        assert_eq!(owner.harmonic().producer_hash(), ContentHash([43; 32]));
    }
}

#[test]
fn weak_measured_wall_work_survives_equal_large_totals_and_remains_signed() {
    let spectrum = spectrum(&segments(), WelchWindow::Rectangular);
    let mut owner = owner(&spectrum);
    owner.impedance = vec![C64::from_re(2e16), C64::ZERO, C64::ZERO, C64::from_re(3e16)];
    owner.exterior_work = owner.impedance.clone();
    owner.wall_work = vec![C64::from_re(2e-8), C64::ZERO, C64::ZERO, C64::from_re(3e-8)];
    let expected = segments()
        .iter()
        .map(|&q| oracle(owner.wall_work_matrix(), q))
        .sum::<f64>()
        / 3.0;
    let mut scratch = [C64::ZERO; 2];
    for sign in [1.0, -1.0] {
        let actual = owner
            .power_welch_bin(&spectrum, BIN, owner.sources(), WORK, &mut scratch, || {
                false
            })
            .unwrap();
        assert_eq!(actual.driving_power_w, actual.exterior_power_w);
        close(actual.absorbed_wall_power_w, sign * expected);
        // The manufactured raw matrices are deliberately imbalanced. Neither
        // the separate weak wall work nor this residual may be repaired away.
        close(actual.balance_residual_w, -sign * expected);
        owner.wall_work.iter_mut().for_each(|z| *z = -*z);
    }
}

#[test]
fn measured_molecular_work_refuses_order_scratch_work_and_cancellation() {
    let spectrum = spectrum(&segments(), WelchWindow::Rectangular);
    let owner = owner(&spectrum);
    let mut scratch = [C64::ZERO; 2];
    for bin in [0, N / 2, BIN + 1] {
        assert!(
            owner
                .power_welch_bin(&spectrum, bin, owner.sources(), WORK, &mut scratch, || {
                    false
                })
                .is_err()
        );
    }
    let reversed = [owner.sources()[1].clone(), owner.sources()[0].clone()];
    assert!(
        owner
            .power_welch_bin(&spectrum, BIN, &reversed, WORK, &mut scratch, || false)
            .is_err()
    );
    assert!(
        owner
            .power_welch_bin(
                &spectrum,
                BIN,
                owner.sources(),
                WORK,
                &mut scratch[..1],
                || false
            )
            .is_err()
    );
    assert!(
        owner
            .power_welch_bin(
                &spectrum,
                BIN,
                owner.sources(),
                HarmonicWorkBudget {
                    maximum_products: 65
                },
                &mut scratch,
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        owner.power_welch_bin(&spectrum, BIN, owner.sources(), WORK, &mut scratch, || {
            polls += 1;
            polls == 6
        }),
        Err(MolecularObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
    let actual = owner
        .power_welch_bin(&spectrum, BIN, owner.sources(), WORK, &mut scratch, || {
            false
        })
        .unwrap();
    close(
        actual.driving_power_w,
        segments()
            .iter()
            .map(|&q| oracle(owner.impedance(), q))
            .sum::<f64>()
            / 3.0,
    );
}

#[test]
fn actual_molecular_boundary_retains_recorded_work_with_the_original_passive_wall() {
    use crate::helmholtz::laplace::complex::ComplexOptions;
    use crate::observation_transfers::{LocalBoundary, LocalWall};
    use crate::panel3d::SpherePanels;
    use fs_phs::WallPin;

    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let panels = surface.areas().len();
    let source = material();
    let medium = Medium {
        density: source.gas().density,
        sound_speed: source.gas().sound_speed,
    };
    let driven: Vec<_> = surface.normals().iter().map(|n| n[2] >= 0.0).collect();
    let area: f64 = surface
        .areas()
        .iter()
        .zip(&driven)
        .filter_map(|(&a, &on)| on.then_some(a))
        .sum();
    let traces = [
        surface
            .normals()
            .iter()
            .zip(&driven)
            .map(|(n, &on)| if on { 1.0 + 0.4 * n[0] } else { 0.0 })
            .collect::<Vec<_>>(),
        driven
            .iter()
            .map(|&on| if on { 1.0 / area } else { 0.0 })
            .collect(),
    ];
    let ids = [
        SourceId::StructuralMode("surface".into()),
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
        driven.iter().map(|&on| usize::from(!on)).collect(),
    )
    .unwrap();
    let model = MolecularBemObservation::new(
        &surface,
        source,
        &ids,
        &traces,
        &[ReceiverSpec {
            id: "physical microphone".into(),
            position_m: [0.0, 0.0, 1.0],
        }],
        0.01,
        [1.0, 10_000.0],
        ComplexOptions {
            max_panels: panels,
            max_bytes: 32 * panels * panels,
            max_work: 512 * (panels * panels) as u64,
        },
    )
    .unwrap()
    .with_local_boundary(wall)
    .unwrap();
    let inputs = segments();
    let spectrum = spectrum(&inputs, WelchWindow::PeriodicHann);
    let prepared = model
        .prepare_response(
            core::f64::consts::TAU * spectrum.frequency_hz(BIN).unwrap(),
            MolecularResponseLimits {
                maximum_retained_entries: 12 + 4 * panels + 16,
                maximum_panel_products: (12 * panels) as u64,
            },
            || false,
        )
        .unwrap();
    let mut scratch = [C64::ZERO; 2];
    let actual = prepared
        .power_welch_bin(&spectrum, BIN, &ids, WORK, &mut scratch, || false)
        .unwrap();
    let reference: Vec<_> = inputs
        .iter()
        .map(|q| prepared.power_w(q, || false).unwrap())
        .collect();
    let scale = (2.0 / 3.0) / 3.0;
    close(
        actual.driving_power_w,
        reference.iter().map(|p| p.driving_power_w).sum::<f64>() * scale,
    );
    close(
        actual.exterior_power_w,
        reference.iter().map(|p| p.exterior_power_w).sum::<f64>() * scale,
    );
    close(
        actual.absorbed_wall_power_w,
        reference
            .iter()
            .map(|p| p.absorbed_wall_power_w)
            .sum::<f64>()
            * scale,
    );
    assert!(actual.absorbed_wall_power_w > 0.0);
    assert_eq!(prepared.source(), model.source());
}
