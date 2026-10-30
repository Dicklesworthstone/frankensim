//! Sealed finite volumes keep their original dynamics without exterior rows.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dg::ApertureDg;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::paired::BladePropagator;
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_couple::dg_port::input::{DgSourceError, DgSourceInput};
use fs_couple::dg_port::{DgBranch, DgPortLoad, DgTerminalChannelInput, DgTerminalChannelKind};
use fs_couple::modal_acoustic_time::{ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingConfig};
use fs_dcontact::Obstacle;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};

const DT: f64 = 1.0 / 96_000.0;
const SEALED: PortId = PortId(7);
const OPEN: PortId = PortId(9);

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn duct(length: f64, elements: usize, outlet: DgBoundary, walls: bool) -> DgDuct {
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length,
            }],
        },
        &gas(),
        DgMesh { elements, order: 2 },
        DgBoundary::Port(PortId(0)),
        outlet,
    )
    .unwrap();
    if walls {
        duct.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
    } else {
        duct
    }
}

fn terminal(legacy: bool) -> PreparedMimoLoad {
    // A positive full matrix has real mutual reaction. The sealed coordinate
    // has zero flow, so only its principal exterior submatrix survives.
    let d = if legacy {
        vec![4e7, 2e6, 3e6, 2e6, 3e7, 5e6, 3e6, 5e6, 2e7]
    } else {
        vec![3e7, 5e6, 5e6, 2e7]
    };
    let n = if legacy { 3 } else { 2 };
    PreparedMimoLoad::new(
        &MimoModel {
            ports: n,
            poles: vec![],
            residues: vec![],
            d,
            e: vec![0.0; n * n],
            reciprocal: true,
        },
        DT,
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 33,
            max_states: 16,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}

fn mixed(legacy: bool, walls: bool) -> DgPortLoad {
    let main = duct(0.05, 3, DgBoundary::Port(PortId(1)), walls)
        .with_side_ports(&[
            DgSidePort {
                id: SEALED,
                downstream_element: 1,
            },
            DgSidePort {
                id: OPEN,
                downstream_element: 2,
            },
        ])
        .unwrap();
    let closed = if legacy {
        DgBoundary::Port(PortId(1))
    } else {
        DgBoundary::Closed
    };
    let ports = if legacy {
        vec![SEALED, OPEN, PortId(1)]
    } else {
        vec![OPEN, PortId(1)]
    };
    let load = DgPortLoad::new_with_terminal_branches(
        main,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        ports,
        terminal(legacy),
        vec![
            DgBranch {
                port: SEALED,
                duct: duct(0.012, 1, closed, walls),
            },
            DgBranch {
                port: OPEN,
                duct: duct(0.018, 1, DgBoundary::Port(PortId(1)), walls),
            },
        ],
    )
    .unwrap();
    if legacy {
        load.with_terminal_channels(
            vec![DgTerminalChannelInput {
                port: SEALED,
                kind: DgTerminalChannelKind::Closed,
            }],
            0,
        )
        .unwrap()
    } else {
        load
    }
}

fn source(load: DgPortLoad, dt: f64) -> ApertureDg {
    let contact = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        0.0,
        2.0,
        "sealed-volume regression: explicit zero contact".into(),
    )
    .unwrap();
    let aperture = DynamicAperture::new_with_propagator(
        DynamicApertureSpec {
            aperture: BernoulliAperture {
                rest_opening_m: 0.0002,
                width_m: 0.001,
                closing_pressure_pa: 10_000.0,
            },
            mass_kg: 1e-5,
            stiffness_n_m: 500.0,
            damping_ratio: 0.0,
            density_kg_m3: load.duct().inlet_gas().density,
            impedance_pa_s_m3: load.impedance_pa_s_m3(),
            time_step_s: dt,
            max_steps: 256,
        },
        ApertureState {
            opening_m: 0.0002,
            opening_velocity_m_s: 0.0,
        },
        contact,
        BladePropagator::ExactZoh,
    )
    .unwrap();
    ApertureDg::new(aperture, load).unwrap()
}

fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 2e-8 * scale.max(1e-25),
        "{a:e} != {b:e} (scale {scale:e})"
    );
}

#[test]
fn sealed_original_volume_matches_zero_flow_termination_with_wall_memory_and_mutual_load() {
    for walls in [false, true] {
        let (a, b) = (mixed(false, walls), mixed(true, walls));
        close(
            a.impedance_pa_s_m3(),
            b.impedance_pa_s_m3(),
            b.impedance_pa_s_m3(),
        );
        assert_eq!(a.terminal_ports(), &[OPEN, PortId(1)]);
        assert_eq!(a.branches()[0].duct().outlet(), DgBoundary::Closed);
        let (mut a, mut b) = (source(a, DT), source(b, DT));
        let mut sealed_peak = 0.0_f64;
        for sample in 0..96 {
            let drive = TubeDrive {
                upstream_pressure_pa: if sample < 24 { 30.0 } else { 0.0 },
                body_flow_m3_s: 0.0,
            };
            a.step(drive).unwrap();
            b.step(drive).unwrap();
            let (a, b) = (a.load(), b.load());
            close(
                a.stored_energy_j(),
                b.stored_energy_j(),
                b.stored_energy_j(),
            );
            for (x, y) in a.state().pressure().iter().zip(b.state().pressure()) {
                close(*x, *y, x.abs() + y.abs());
            }
            for (x, y) in a.branches().iter().zip(b.branches()) {
                close(
                    x.state().stored_energy_j(),
                    y.state().stored_energy_j(),
                    y.state().stored_energy_j(),
                );
                for (p, q) in x.state().pressure().iter().zip(y.state().pressure()) {
                    close(*p, *q, p.abs() + q.abs());
                }
            }
            let sealed = &a.last_branch_frames()[0];
            sealed_peak = sealed_peak.max(a.branches()[0].state().stored_energy_j());
            assert_eq!(sealed.duct.port_outputs.len(), 1);
            assert!(
                !a.last_terminal_frames()
                    .iter()
                    .any(|frame| frame.port == SEALED)
            );
            close(
                sealed.bore_work_j,
                sealed.branch_work_j,
                sealed.bore_work_j.abs() + sealed.branch_work_j.abs(),
            );
            for (actual, expected) in a
                .last_terminal_frames()
                .iter()
                .zip(&b.last_terminal_frames()[1..])
            {
                assert_eq!(actual.port, expected.port);
                close(
                    actual.outward_flow_m3_s,
                    expected.outward_flow_m3_s,
                    actual.outward_flow_m3_s.abs() + expected.outward_flow_m3_s.abs(),
                );
            }
        }
        assert!(
            sealed_peak > 1e-16,
            "a sealed chimney must store actual acoustic energy"
        );
    }
}

fn sealed_load(length: f64) -> DgPortLoad {
    let main = duct(0.05, 2, DgBoundary::Closed, false)
        .with_side_ports(&[DgSidePort {
            id: SEALED,
            downstream_element: 1,
        }])
        .unwrap();
    DgPortLoad::new_with_closed_branches(
        main,
        DT,
        1e-11,
        64,
        vec![DgBranch {
            port: SEALED,
            duct: duct(length, 1, DgBoundary::Closed, false),
        }],
    )
    .unwrap()
}

#[test]
fn sealed_volume_changes_source_feedback_and_continuous_resonance_without_an_exterior() {
    let mut driven = source(sealed_load(0.012), DT);
    for _ in 0..32 {
        driven
            .step(TubeDrive {
                upstream_pressure_pa: 20.0,
                body_flow_m3_s: 0.0,
            })
            .unwrap();
    }
    assert!(driven.load().terminal().is_none());
    assert!(driven.load().terminal_ports().is_empty());
    assert!(driven.load().branches()[0].state().stored_energy_j() > 0.0);
    assert!(
        driven
            .load()
            .port_response()
            .unwrap()
            .history_pressure_pa
            .abs()
            > 0.0
    );

    let gate = CancelGate::new_clock_free();
    let mechanics = CoupledModalSystem::new(
        vec![
            ModalAcousticTimeModel::try_free_mass(
                96_000,
                0.1,
                0.0,
                0.0,
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
        ],
        vec![],
        ModalCouplingConfig {
            max_modes: 4,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1000.0,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-11,
            energy_relative_tolerance: 1e-9,
        },
        &gate,
    )
    .unwrap();
    let response = |length| {
        mechanics
            .harmonic_response_with_load(
                6000.0,
                &[C64::from_re(0.01)],
                &sealed_load(length),
                &[1e-5],
                10_000_000,
                &gate,
            )
            .unwrap()
    };
    let (short, long) = (response(0.012), response(0.024));
    assert!(short.terminal.is_none());
    assert_eq!(short.branches.len(), 1);
    assert_eq!(short.branches[0].duct.ports.len(), 1);
    assert!(
        (short.inlet_pressure_pa - long.inlet_pressure_pa).abs()
            > 1e-4 * short.inlet_pressure_pa.abs()
    );
    let branch = &short.branches[0];
    close(
        branch.bore_power_w,
        branch.branch_power_w,
        branch.bore_power_w.abs() + branch.branch_power_w.abs() + 1e-15,
    );
}

fn authored_chimneys() -> String {
    let mut text = String::from(
        "dg_tube 0.05 0.001 baffled 10 1 8 50 2000 0.01 1e-11 1024\ndg_side_openings 9\n",
    );
    for index in 1..=9 {
        text.push_str(&format!(
            "dg_side_chimney {} {index} {} 0 0.001 0.012 1 1\ndg_chimney_gas 293.15 101325 0\n",
            index + 1,
            0.01 * index as f64
        ));
    }
    text.push_str("dg_terminal_channels 8 0\n");
    for port in 2..=9 {
        text.push_str(&format!("dg_terminal_closed {port}\n"));
    }
    text.push_str("dg_radiation computed-baffled 20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257\ndg_exterior 293.15 101325 0\n");
    text
}

fn read(text: &str) -> Result<DgSourceInput, DgSourceError> {
    let mut lines = text.lines().enumerate();
    let (line, header) = lines.next().unwrap();
    let mut fields = header.split_ascii_whitespace();
    let kind = fields.next().unwrap();
    DgSourceInput::read(
        kind,
        line + 1,
        &mut fields,
        &mut || {
            let (line, row) = lines.next().ok_or(DgSourceError::Input {
                line: 0,
                what: "missing row",
            })?;
            Ok((line + 1, row.split_ascii_whitespace()))
        },
        &gas(),
        None,
        &mut Vec::new(),
        &CancelGate::new_clock_free(),
    )
}

#[test]
fn authored_full_chimney_inventory_bakes_and_observes_only_real_exterior_openings() {
    let gate = CancelGate::new_clock_free();
    let load = read(&authored_chimneys())
        .unwrap()
        .prepare(&gas(), 1.0 / 8000.0, &gate)
        .unwrap();
    assert_eq!(load.branches().len(), 9);
    assert_eq!(load.terminal_ports(), &[PortId(1), PortId(10)]);
    assert_eq!(load.radiation_record().unwrap().apertures().len(), 2);
    assert!(load.terminal_channels().is_empty());
    assert_eq!(
        load.branches()
            .iter()
            .filter(|b| b.duct().outlet() == DgBoundary::Closed)
            .count(),
        8
    );
    let mut observer = load
        .baffled_receiver(
            8000,
            CircularOutletReceiver {
                position_m: [0.1, 0.05, 0.2],
                radial_rings: 2,
                angular_points: 8,
                maximum_frequency_hz: 1000.0,
            },
        )
        .unwrap();
    let mut model = source(load, 1.0 / 8000.0);
    for _ in 0..32 {
        model
            .step(TubeDrive {
                upstream_pressure_pa: 1.0,
                body_flow_m3_s: 0.0,
            })
            .unwrap();
        let flows: Vec<_> = model
            .load()
            .last_terminal_frames()
            .iter()
            .map(|f| f.outward_flow_m3_s)
            .collect();
        assert_eq!(flows.len(), 2);
        assert!(observer.step_trace(&flows).unwrap().is_finite());
    }
    assert!(
        model.load().branches()[..8]
            .iter()
            .all(|b| b.state().step_count() == 32)
    );
    assert!(
        model.load().branches()[..8]
            .iter()
            .any(|b| b.state().stored_energy_j() > 0.0)
    );
}

#[test]
fn actual_exterior_and_complete_internal_state_budgets_still_refuse() {
    let text = authored_chimneys();
    // Keep only one seal: nine actual openings still exceed the eight-row owner.
    let mut radiating = text.replace("dg_terminal_channels 8 0", "dg_terminal_channels 1 0");
    for port in 3..=9 {
        radiating = radiating.replace(&format!("dg_terminal_closed {port}\n"), "");
    }
    assert!(
        matches!(read(&radiating), Err(DgSourceError::Input { what, .. })
        if what == "unsealed DG openings exceed the exterior radiation aperture budget")
    );
    assert!(read(&text.replace("1e-11 1024", "1e-11 400")).is_err());
    assert!(read(&text.replace("dg_side_openings 9", "dg_side_openings 33")).is_err());
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(
        read(&text)
            .unwrap()
            .prepare(&gas(), 1.0 / 8000.0, &gate)
            .is_err()
    );
}
