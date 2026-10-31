//! Analytic modal authority, separate spatial/time ladders, and actual runtime
//! work witnesses for the explicitly qualified lossless higher-mode line.
use core::f64::consts::{PI, TAU};
use fs_couple::driving_point::LineLadder;
use fs_couple::mm_line::evanescent::{
    EvanescentConfig, EvanescentError, MatchedEvanescentPorts, UniformModalLine,
};
use fs_couple::mm_line::{MmLineBank, MmLineConfig, MmLoad};
use fs_duct::modal::{j1_roots, modal_characteristic_impedances};
use fs_duct::{Duct, LossModel, Segment, TabulatedLoad, TabulatedMouth, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn cutoff(gas: &GasState, radius: f64) -> f64 {
    gas.sound_speed * j1_roots(1).unwrap()[0] / radius
}
fn config(cells: usize, band: f64, error: f64) -> EvanescentConfig {
    EvanescentConfig {
        cells_per_mode: cells,
        band_max_rad_s: band,
        max_reflection_error: error,
    }
}
fn reflection(z: C64, z0: f64) -> C64 {
    (z - C64::from_re(z0)) / (z + C64::from_re(z0))
}
fn source(gas: &GasState, radius: f64, omega: f64, mode: usize) -> C64 {
    // Native frequency authority is exp(-i wt); the causal realizer is exp(+i wt).
    modal_characteristic_impedances(gas, radius, omega, LossModel::Lossless, mode + 1).unwrap()
        [mode]
        .conj()
}

#[test]
fn three_quadrature_rungs_converge_to_native_evanescent_impedance() {
    let gas = air();
    let radius = 0.01;
    let wc = cutoff(&gas, radius);
    let dt = 0.001 / wc;
    let band = 0.9 * wc;
    let mut relative = Vec::new();
    for n in [2, 4, 8] {
        let ports = MatchedEvanescentPorts::new(&gas, radius, 2, dt, config(n, band, 1.0)).unwrap();
        let cells = ports.cell_spec(0).unwrap();
        assert_eq!(cells.cells().len(), n);
        assert!(
            cells
                .cells()
                .iter()
                .all(|cell| cell.resistance_pa_s_m3 == f64::INFINITY)
        );
        assert!((cells.static_inertance_pa_s2_m3() / (ports.zc0() / wc) - 1.0).abs() < 4e-15);
        let reference = source(&gas, radius, band, 1);
        let actual = cells.impedance(C64::new(0.0, band));
        relative.push((actual - reference).abs() / reference.abs());
        assert_eq!(actual.re, 0.0);
        let receipt = &ports.receipts()[0];
        for k in 1..=64 {
            let omega = band * f64::from(k) / 64.0;
            let warped = 2.0 * (0.5 * omega * dt).tan() / dt;
            let expected = reflection(source(&gas, radius, omega, 1), ports.zc0());
            let response = reflection(cells.impedance(C64::new(0.0, warped)), ports.zc0());
            assert!((response - expected).abs() <= receipt.reflection_error_bound + 3e-12);
            assert!((response.abs() - 1.0).abs() < 3e-14);
        }
    }
    // Distinct, live spatial errors; doubling N must improve the actual owner.
    assert!(relative[0] > 0.04 && relative[0] < 0.05);
    assert!(relative[1] > 0.001 && relative[1] < 0.002);
    assert!(relative[2] > 6e-7 && relative[2] < 7e-7);
}

#[test]
fn three_clock_rungs_bound_complex_reflection_error_separately() {
    let gas = air();
    let radius = 0.01;
    let wc = cutoff(&gas, radius);
    let band = 0.7 * wc;
    let mut errors = Vec::new();
    for step in [0.6, 0.3, 0.15] {
        let dt = step / wc;
        let ports =
            MatchedEvanescentPorts::new(&gas, radius, 3, dt, config(16, band, 0.1)).unwrap();
        for (i, receipt) in ports.receipts().iter().enumerate() {
            let cells = ports.cell_spec(i).unwrap();
            for k in 1..=64 {
                let omega = band * f64::from(k) / 64.0;
                let warped = 2.0 * (0.5 * omega * dt).tan() / dt;
                let actual = reflection(cells.impedance(C64::new(0.0, warped)), ports.zc0());
                let expected = reflection(source(&gas, radius, omega, i + 1), ports.zc0());
                assert!((actual - expected).abs() <= receipt.reflection_error_bound + 3e-12);
                if i == 0 && k == 64 {
                    errors.push((actual - expected).abs());
                }
            }
            assert!(receipt.quadrature_reflection_error_bound < 1e-10);
        }
    }
    assert!(errors[0] > 0.01);
    for pair in errors.windows(2) {
        assert!(
            pair[1] > 0.2 * pair[0] && pair[1] < 0.3 * pair[0],
            "{errors:?}"
        );
    }
}

#[test]
fn actual_sampled_harmonic_matches_the_physical_cell_transfer() {
    let gas = air();
    let radius = 0.01;
    let wc = cutoff(&gas, radius);
    let dt = 0.2 / wc;
    let mut ports =
        MatchedEvanescentPorts::new(&gas, radius, 2, dt, config(4, 0.4 * wc, 0.01)).unwrap();
    let samples = 4096;
    let theta = TAU * 37.0 / samples as f64;
    let warped = 2.0 * (theta / 2.0).tan() / dt;
    let expected = reflection(
        ports.cell_spec(0).unwrap().impedance(C64::new(0.0, warped)),
        ports.zc0(),
    );
    let mut measured = C64::ZERO;
    for n in 0..9 * samples {
        let phase = theta * n as f64;
        let a = phase.cos();
        let affine = ports.boundary_response(0).unwrap();
        let b = ports.push(&[a]).unwrap()[0];
        assert!((b - affine.incoming(a).unwrap()).abs() < 2e-13);
        if n >= 8 * samples {
            measured = measured + C64::new(phase.cos(), -phase.sin()).scale(b);
        }
    }
    measured = measured.scale(2.0 / samples as f64);
    assert!(
        (measured - expected).abs() < 2e-10,
        "{measured:?} != {expected:?}"
    );
}

#[test]
fn all_higher_ports_obey_unmodified_storage_and_wave_work_for_arbitrary_inputs() {
    let gas = air();
    let radius = 0.01;
    let wc = cutoff(&gas, radius);
    let dt = 0.2 / wc;
    let mut ports =
        MatchedEvanescentPorts::new(&gas, radius, 8, dt, config(8, 0.6 * wc, 0.02)).unwrap();
    let z0 = ports.zc0();
    let mut wave_work = 0.0;
    let mut absolute_work = 0.0;
    let mut released = false;
    for n in 0..4096 {
        let a: [f64; 7] = core::array::from_fn(|j| {
            if n < 200 {
                100.0 * ((n * (j + 3) % 19) as f64 - 9.0)
            } else {
                0.0
            }
        });
        let b = ports.push(&a).unwrap().to_vec();
        let mut independent_storage = 0.0;
        for (j, frame) in ports.frames().iter().enumerate() {
            let work = dt * (a[j] * a[j] - b[j] * b[j]) / z0;
            wave_work += work;
            absolute_work += work.abs();
            released |= n >= 200 && work < 0.0;
            assert_eq!(frame.dissipated_energy_j, 0.0);
            let scale =
                (frame.stored_energy_j + frame.supplied_work_j.abs()).max(f64::MIN_POSITIVE);
            assert!(frame.balance_residual_j().abs() < 4e-12 * scale);
            assert!(
                (frame.pressure_pa - (a[j] + b[j])).abs() < 3e-12 * (1.0 + a[j].abs() + b[j].abs())
            );
            assert!((frame.flow_m3_s - (a[j] - b[j]) / z0).abs() < 3e-15);
            for (k, cell) in ports.cell_spec(j).unwrap().cells().iter().enumerate() {
                independent_storage +=
                    0.5 * cell.inertance_pa_s2_m3 * frame.inductor_flows_m3_s()[k].powi(2)
                        + 0.5 * cell.compliance_m3_pa * frame.pressures_pa()[k].powi(2);
            }
        }
        assert!(
            (ports.stored_energy_j() - independent_storage).abs()
                < 3e-14 * independent_storage.max(f64::MIN_POSITIVE)
        );
        assert!(
            (wave_work - independent_storage).abs() < 5e-12 * absolute_work.max(f64::MIN_POSITIVE)
        );
    }
    assert!(
        released,
        "stored reactive energy must be returned without clipping"
    );
}

#[test]
fn a_late_mode_refusal_preserves_every_accepted_state_and_observation() {
    let gas = air();
    let wc = cutoff(&gas, 0.01);
    let make =
        || MatchedEvanescentPorts::new(&gas, 0.01, 4, 0.2 / wc, config(4, 0.5 * wc, 0.02)).unwrap();
    let (mut ports, mut control) = (make(), make());
    for a in [[1.0, 3.0, -2.0], [-4.0, 0.5, 9.0]] {
        ports.push(&a).unwrap();
        control.push(&a).unwrap();
    }
    let old_frames = ports.frames().to_vec();
    let old_incoming = ports.incoming().to_vec();
    let old_energy = ports.stored_energy_j();
    assert!(ports.push(&[1.0, f64::MAX, 0.0]).is_err());
    assert!(ports.push(&[f64::NAN, 0.0, 0.0]).is_err());
    assert_eq!(ports.frames(), old_frames);
    assert_eq!(ports.incoming(), old_incoming);
    assert_eq!(ports.stored_energy_j(), old_energy);
    assert_eq!(
        ports.push(&[0.2, 0.7, -0.1]).unwrap(),
        control.push(&[0.2, 0.7, -0.1]).unwrap()
    );
}

fn uniform(rate: u32) -> (Duct, GasState, TabulatedLoad, MmLineConfig) {
    let gas = air();
    let radius = 0.002;
    let length = 8.0 * gas.sound_speed / f64::from(rate);
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let z0 = gas.density * gas.sound_speed / (PI * radius * radius);
    let dt = 1.0 / f64::from(rate);
    let table = TabulatedLoad::try_new(
        vec![TAU / (LineLadder::default().max_taps as f64 * dt), PI / dt],
        vec![C64::from_re(0.7 * z0); 2],
        "constant-resistance analytic control",
        0.0,
    )
    .unwrap()
    .with_mouth(TabulatedMouth {
        axial_position_m: length,
        radius_m: radius,
        gas,
        dc_impedance: 0.7 * z0,
    })
    .unwrap();
    (
        duct,
        gas,
        table,
        MmLineConfig {
            sample_rate_hz: rate,
            n_modes: 2,
            extra_slices: 1,
            loss: LossModel::Lossless,
        },
    )
}

#[test]
fn prior_uniform_fir_refusal_now_has_a_qualified_causal_runtime_and_unchanged_plane() {
    for rate in [48_000, 96_000] {
        let (duct, gas, table, cfg) = uniform(rate);
        let mut line = UniformModalLine::new(
            &duct,
            &gas,
            MmLoad::Tabulated(&table),
            &cfg,
            config(4, TAU * 6000.0, 0.01),
        )
        .unwrap();
        let mut plane = MmLineBank::new(
            core::slice::from_ref(&duct),
            &["scalar control"],
            &gas,
            MmLoad::Tabulated(&table),
            &MmLineConfig { n_modes: 1, ..cfg },
        )
        .unwrap();
        assert!((line.plane_receipt().round_trip_samples - 16.0).abs() < 1e-12);
        assert!(line.higher_ports().receipts()[0].reflection_error_bound < 0.01);
        for n in 0..512 {
            let a = [
                if n < 32 {
                    1.0 + 0.15 * (0.21 * f64::from(n)).cos()
                } else {
                    0.0
                },
                if n < 16 {
                    20.0 * (0.19 * f64::from(n)).cos()
                } else {
                    0.0
                },
            ];
            let response = [
                line.boundary_response(0).unwrap(),
                line.boundary_response(1).unwrap(),
            ];
            let b = line.push(&a).unwrap();
            assert_eq!(
                b[0],
                plane.push(a[0]).unwrap(),
                "original scalar arithmetic must remain identical"
            );
            for j in 0..2 {
                assert!((b[j] - response[j].incoming(a[j]).unwrap()).abs() < 2e-12);
            }
            if n == 0 {
                assert!(
                    b[1].abs() > 1.0,
                    "reactive feedthrough cannot acquire an extra sample delay"
                );
            }
        }
        let old = line.incoming().to_vec();
        let pending = [
            line.boundary_response(0).unwrap(),
            line.boundary_response(1).unwrap(),
        ];
        assert!(line.push(&[1.0, f64::MAX]).is_err());
        assert_eq!(line.incoming(), old);
        assert_eq!(
            [
                line.boundary_response(0).unwrap(),
                line.boundary_response(1).unwrap()
            ],
            pending
        );
    }
}

#[test]
fn domains_and_error_budgets_refuse_instead_of_changing_the_source() {
    let gas = air();
    let wc = cutoff(&gas, 0.01);
    assert!(matches!(
        MatchedEvanescentPorts::new(&gas, 0.01, 2, 1.0 / wc, config(8, 0.99 * wc, 1.0)),
        Err(EvanescentError::Invalid(_))
    ));
    assert!(matches!(
        MatchedEvanescentPorts::new(&gas, 0.01, 2, 0.6 / wc, config(2, 0.7 * wc, 1e-9)),
        Err(EvanescentError::Accuracy { mode: 1, .. })
    ));
    for n in [0, 17] {
        assert!(
            MatchedEvanescentPorts::new(&gas, 0.01, 2, 0.2 / wc, config(n, 0.6 * wc, 1.0)).is_err()
        );
    }
    let (duct, gas, table, cfg) = uniform(48_000);
    let approximation = config(4, TAU * 6000.0, 0.01);
    for load in [
        MmLoad::Analytic(Termination::Closed),
        MmLoad::Analytic(Termination::IdealOpen),
        MmLoad::FlangedModal(&gas),
        MmLoad::Matched,
    ] {
        assert!(matches!(
            UniformModalLine::new(&duct, &gas, load, &cfg, approximation),
            Err(EvanescentError::Invalid(_))
        ));
    }
    assert!(matches!(
        UniformModalLine::new(
            &duct,
            &gas,
            MmLoad::Tabulated(&table),
            &MmLineConfig {
                loss: LossModel::Bessel,
                ..cfg
            },
            approximation
        ),
        Err(EvanescentError::Invalid(_))
    ));
    let nonuniform = Duct {
        segments: vec![
            duct.segments[0].clone(),
            Segment::Cylinder {
                radius: 0.003,
                length: 0.1,
            },
        ],
    };
    assert!(matches!(
        UniformModalLine::new(
            &nonuniform,
            &gas,
            MmLoad::Tabulated(&table),
            &cfg,
            approximation
        ),
        Err(EvanescentError::Invalid(_))
    ));
}
