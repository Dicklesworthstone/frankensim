//! G0 admitted schema/accepted-account laws; G3 independently flipped work.
//! Existing synthetic material fixtures remain source-backed, not recordings.
use fs_couple::bernoulli_aperture::{
    BernoulliAperture,
    dynamic::{ApertureDrive, ApertureState, DynamicAperture, DynamicApertureSpec},
    modal::{ApertureMode, ModalAperture, ModalApertureSpec},
    paired::{BladeState, PairedApertureLimits},
};
use fs_couple::cavity_modal::shared::{CavitySurface, SharedCavity};
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::cavity_phs::piston::double_acting::*;
use fs_couple::cavity_phs::piston::*;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalColumnConnection, ModalCouplingConfig,
};
use fs_couple::*;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_math::c64::C64;
use fs_qty::Energy;

// Reuse the existing public-API synthetic gas fixture and its regressions.
// No production constructor gains a preset or test-only material fallback.
#[path = "../src/cavity_phs/piston/tests.rs"]
mod piston_fixture;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn budget() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 4,
        max_setup_terms: 1024,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 1.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-8,
    }
}
fn model(w: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        48000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: w,
            damping_ratio: 0.01,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}
fn coupled() -> CoupledModalSystem {
    CoupledModalSystem::new_with_columns(
        vec![model(400.), model(600.), model(800.)],
        vec![ModalColumnConnection {
            column: vec![1., -2., 0.5],
            stiffness_n_m: 10.,
            damping_n_s_m: 0.02,
            rest_extension_m: 0.,
        }],
        budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn shared() -> SharedCavity {
    let source = HelmholtzModalCavity::new(
        fs_scenario::HelmholtzCavity {
            volume_m3: 0.01,
            neck_radius_m: 0.01,
            neck_length_m: 0.02,
            inner_end: fs_scenario::NeckEnd::Flanged,
            outer_end: fs_scenario::NeckEnd::Flanged,
        },
        air(),
        CavityNeckLoss::Resistance(2000.),
    )
    .unwrap();
    SharedCavity::new(
        source,
        vec![
            CavitySurface {
                model: model(500.),
                outward_volume_shapes: vec![0.01],
            },
            CavitySurface {
                model: model(900.),
                outward_volume_shapes: vec![-0.02],
            },
        ],
        48000,
        ModalAcousticState::default(),
        ModalAcousticTimeBudget::audible_reference(),
        budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn obstacle() -> Obstacle {
    Obstacle::new(
        vec![-1.],
        1,
        1,
        vec![0.],
        vec![1.],
        1e8,
        2.,
        "synthetic adoption obstacle".into(),
    )
    .unwrap()
}
fn dynamic() -> DynamicAperture {
    DynamicAperture::new(
        DynamicApertureSpec {
            aperture: BernoulliAperture {
                rest_opening_m: 0.001,
                width_m: 0.01,
                closing_pressure_pa: 10000.,
            },
            mass_kg: 0.001,
            stiffness_n_m: 1000.,
            damping_ratio: 0.1,
            density_kg_m3: 1.2,
            impedance_pa_s_m3: 1e6,
            time_step_s: 1. / 192000.,
            max_steps: 32,
        },
        ApertureState {
            opening_m: 0.001,
            opening_velocity_m_s: 0.,
        },
        obstacle(),
    )
    .unwrap()
}
fn modal() -> ModalAperture {
    ModalAperture::new(
        ModalApertureSpec {
            rest_opening_m: 0.001,
            width_m: 0.01,
            modes: vec![ApertureMode {
                mass_kg: 0.001,
                stiffness_n_m: 1000.,
                damping_n_s_m: 0.01,
                pressure_area_m2: 1e-4,
                gap_projection: 1.,
                force_projection: [1., 0.5],
            }],
            density_kg_m3: 1.2,
            impedance_pa_s_m3: 1e6,
            time_step_s: 1. / 192000.,
            max_steps: 32,
            limits: PairedApertureLimits {
                maximum_abs_displacement_m: 0.01,
                maximum_abs_velocity_m_s: 100.,
                maximum_energy_j: 1.,
                max_evaluations: 256,
                velocity_tolerance_m_s: 1e-13,
                relative_tolerance: 1e-11,
            },
        },
        vec![BladeState::default()],
        obstacle(),
    )
    .unwrap()
}
fn reed() -> fs_couple::render::ReedBoreVoice {
    fs_couple::render::ReedBoreVoice::new(
        &fs_duct::Duct {
            segments: vec![fs_duct::Segment::Cylinder {
                radius: 0.0022,
                length: 0.5,
            }],
        },
        &air(),
        fs_scenario::BeatingReed {
            rest_opening_m: 4e-4,
            width_m: 0.013,
            closing_pressure_pa: 6000.,
            blowing_pressure_pa: 2800.,
            attack_s: 0.008,
            mass_kg: 0.,
            stiffness_n_m: 0.,
            damping_ratio: 0.35,
        },
        fs_duct::Termination::UnflangedOpen,
        fs_couple::thin_plate::PlateBank::default(),
        1.,
        48000,
        4800,
        None,
    )
    .unwrap()
}
fn single() -> GasPiston {
    GasPiston::new(piston_fixture::gas(0.75, 0.001), 1., 0.01, 0.).unwrap()
}
fn double() -> DoubleActingGasPiston {
    DoubleActingGasPiston::new(single(), piston_fixture::gas(0.25, 0.001), 0.01).unwrap()
}
fn double_controls() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: piston_fixture::controls(),
        opposing_gas: piston_fixture::controls().gas,
        volume_tolerance_m3: 1e-18,
    }
}
fn drive() -> ApertureDrive {
    ApertureDrive {
        upstream_pressure_pa: 20.,
        incoming_pressure_pa: 1.,
        body_flow_m3_s: 1e-7,
    }
}

fn schemas(module: &dyn PortedModule) {
    assert!(!module.port_schemas().is_empty());
    assert_eq!(
        module.port_coordinate_charts().len(),
        module.port_schemas().len()
    );
    for s in module.port_schemas() {
        assert_eq!(
            PortSchema::try_new(
                s.id().clone(),
                s.kind(),
                s.effort_dimensions(),
                s.flow_dimensions(),
                s.shape(),
                s.coordinates().clone(),
                s.power_pairing(),
                s.timestamp().clone(),
                s.conservation_roles().iter().copied()
            )
            .unwrap(),
            *s
        );
        assert_eq!(
            s.effort_dimensions()
                .checked_plus(s.flow_dimensions())
                .unwrap(),
            fs_qty::Power::DIMS
        );
        for role in s.kind().required_conservation_roles() {
            assert!(s.conservation_roles().contains(role));
        }
    }
    assert!(module.last_window_transfers().is_none());
    assert!(module.last_window_transfer_error().is_none());
    assert_eq!(
        module.audit_last_window().unwrap(),
        ModuleWindowAudit::NotStepped
    );
}
fn receipt(name: &str, module: &dyn PortedModule, available: bool) {
    let transfers = module.last_window_transfers().expect("accepted transfers");
    assert_eq!(transfers.len(), module.port_schemas().len());
    assert!(module.last_window_duration_s().unwrap() > 0.);
    assert!(module.last_window_transfer_error().is_none());
    let sha = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap();
    assert!(sha.status.success());
    let sha = String::from_utf8(sha.stdout).unwrap();
    let sha = sha.trim();
    match module.audit_last_window().unwrap() {
        ModuleWindowAudit::Available(report) => {
            assert!(available);
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"port-adoption\",\"module\":\"{name}\",\"ports\":{},\"audit_violations\":{},\"energy_residual_j\":{},\"tolerance_j\":{},\"verdict\":\"{}\",\"git_sha\":\"{sha}\"}}",
                transfers.len(),
                report.violations().len(),
                report.residuals().first_law().value(),
                report.spec().tolerances().energy().value(),
                if report.violations().is_empty() {
                    "pass"
                } else {
                    "fail"
                }
            );
            assert!(
                report.violations().is_empty(),
                "{name}: {:?}",
                report.violations()
            );
        }
        ModuleWindowAudit::Unavailable { owner, reason } => {
            assert!(!available);
            assert!(!owner.as_str().is_empty());
            assert!(reason.contains("energy tolerance"));
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"port-adoption\",\"module\":\"{name}\",\"ports\":{},\"audit_violations\":null,\"energy_residual_j\":null,\"tolerance_j\":null,\"verdict\":\"unavailable\",\"git_sha\":\"{sha}\",\"reason\":\"{reason}\"}}",
                transfers.len()
            );
        }
        other => panic!("accepted module has no account: {other:?}"),
    }
}

#[test]
fn five_modules_expose_admitted_schemas() {
    let (c, s, p, d, a, m, r) = (
        coupled(),
        shared(),
        single(),
        double(),
        dynamic(),
        modal(),
        reed(),
    );
    for module in [&c as &dyn PortedModule, &s, &p, &d, &a, &m, &r] {
        schemas(module);
    }
    assert_eq!(c.port_schemas().len(), 1);
    assert_eq!(s.port_schemas().len(), 3);
    assert!(matches!(
        s.port_coordinate_charts()[0],
        PortCoordinateChart::MassNormalized {
            reference_mass_kg: 1.,
            volume_projection_m2_per_sqrt_kg: Some(0.01)
        }
    ));
    assert_eq!(
        p.port_schemas()[1].kind(),
        PortKind::ThermalTemperatureEntropy
    );
    assert_eq!(
        d.port_schemas()
            .iter()
            .filter(|p| p.kind() == PortKind::ThermalTemperatureEntropy)
            .count(),
        2
    );
    assert_eq!(m.port_schemas().len(), 4);
}

#[test]
fn window_audit_closes_per_module() {
    let mut c = coupled();
    c.step(&[10., -20., 5.]).unwrap();
    receipt("CoupledModalSystem", &c, true);
    let mut s = shared();
    s.step(&[10., -20.], 2., &CancelGate::new_clock_free())
        .unwrap();
    receipt("SharedCavity", &s, true);
    let mut p = single();
    p.step(
        piston_fixture::drive(20.),
        piston_fixture::controls(),
        || false,
    )
    .unwrap();
    receipt("GasPiston", &p, true);
    let mut d = double();
    d.step(
        DoublePistonInput {
            // The fixture's resolved clock: about 120 rad/s refuses at 0.0005 s.
            duration_s: piston_fixture::drive(0.).duration_s,
            heat_j: [20., 0.],
            external_force_n: 0.,
        },
        double_controls(),
        || false,
    )
    .unwrap();
    receipt("DoubleActingGasPiston", &d, true);
    let mut a = dynamic();
    let f = a.step_with_force(drive(), 0.001).unwrap();
    receipt("DynamicAperture", &a, false);
    let work: f64 = a
        .last_window_transfers()
        .unwrap()
        .iter()
        .map(|t| t.energy_into().value())
        .sum();
    assert!(
        (work - f.pressure_work_j - f.mechanical_work_j).abs()
            <= 16. * f64::EPSILON * work.abs().max(f64::MIN_POSITIVE)
    );
    let mut m = modal();
    m.step(drive(), [0.001, -0.0002]).unwrap();
    receipt("ModalAperture", &m, true);
    let mut r = reed();
    r.step_block(&mut [0.; 2]).unwrap();
    receipt("ReedBoreVoice", &r, false);
}

#[test]
fn sign_flipped_adapter_is_caught() {
    let mut s = shared();
    s.step(&[10., -20.], 2., &CancelGate::new_clock_free())
        .unwrap();
    let ModuleWindowAudit::Available(report) = s.audit_last_window().unwrap() else {
        panic!("audit unavailable")
    };
    assert!(report.violations().is_empty());
    let mut flipped = false;
    let contributions = report
        .contributions()
        .iter()
        .map(|row| {
            if !flipped && row.manifest_entry().contribution_id().as_str() == "port/0" {
                flipped = true;
                let t = row.transfer();
                let bad = IntegratedWindowTransfer::try_new(
                    Energy::new(-t.energy_into().value()),
                    t.mass_into(),
                    t.element_amounts_into().iter().cloned(),
                    t.charge_into(),
                    t.entropy_into().clone(),
                    t.entropy_convention().clone(),
                )
                .unwrap();
                WindowBoundaryContribution::try_new(
                    report.spec(),
                    row.manifest_entry().contribution_id(),
                    bad,
                    row.integration_evidence().clone(),
                )
                .unwrap()
            } else {
                row.clone()
            }
        })
        .collect::<Vec<_>>();
    assert!(flipped);
    let bad = WindowAuditReport::audit(
        report.spec().clone(),
        report.initial().clone(),
        report.final_snapshot().clone(),
        contributions,
    )
    .unwrap();
    assert!(!bad.violations().is_empty());
    assert!(bad.residuals().first_law().value().abs() > bad.spec().tolerances().energy().value());
}

#[test]
fn n_port_junction_two_port_instance_equals_conservative_junction() {
    let c = coupled();
    let base = &c.junction_schemas()[0];
    let ports = base.ports()[..2].to_vec();
    let id = StableId::new("two-port").unwrap();
    assert_eq!(
        JunctionSchema::try_new(
            id.clone(),
            ports.clone(),
            JunctionLaw::CommonEffort {
                flow_signs: vec![1, 1]
            }
        )
        .map(|_| ()),
        ConservativeJunction::new(id, ports[0].clone(), ports[1].clone()).map(|_| ())
    );
    let id = ports[0].id().clone();
    assert_eq!(
        JunctionSchema::try_new(
            id.clone(),
            ports.clone(),
            JunctionLaw::CommonEffort {
                flow_signs: vec![1, 1]
            }
        )
        .map(|_| ()),
        ConservativeJunction::new(id, ports[0].clone(), ports[1].clone()).map(|_| ())
    );
}
#[test]
fn three_component_modal_network_uses_one_n_port_junction() {
    let c = coupled();
    assert_eq!(c.junction_schemas().len(), 1);
    assert_eq!(c.junction_schemas()[0].ports().len(), 4); // Three bodies plus connection storage.
    let s = shared();
    assert_eq!(s.junction_schemas().len(), 1);
    assert_eq!(s.junction_schemas()[0].ports().len(), 4); // Two surfaces, neck, compression storage.
}
#[test]
fn double_acting_piston_faces_share_one_wall_junction() {
    let mut d = double();
    assert_eq!(d.junction_schemas().len(), 1);
    assert_eq!(d.junction_schemas()[0].ports().len(), 3);
    assert!(matches!(
        d.junction_schemas()[0].law(),
        JunctionLaw::CommonFlow { .. }
    ));
    d.step(
        DoublePistonInput {
            duration_s: piston_fixture::drive(0.).duration_s,
            heat_j: [20., 0.],
            external_force_n: 0.,
        },
        double_controls(),
        || false,
    )
    .unwrap();
    let t = d.last_window_transfers().unwrap();
    let net = t[..3].iter().map(|t| t.energy_into().value()).sum::<f64>();
    assert!(net.abs() <= double_controls().piston.energy_tolerance_j);
    assert!(t[0].energy_into().value() > 0.);
    assert!(t[1].energy_into().value() < 0.);
}
#[test]
fn junction_refuses_non_conjugate_or_non_skew_incidence() {
    let c = coupled();
    let mut ports = c.junction_schemas()[0].ports()[..2].to_vec();
    let fluid = reed();
    ports[1] = fluid.port_schemas()[0].clone();
    assert!(
        JunctionSchema::try_new(
            StableId::new("bad-kind").unwrap(),
            ports,
            JunctionLaw::CommonEffort {
                flow_signs: vec![1, 1]
            }
        )
        .is_err()
    );
    let ports = c.junction_schemas()[0].ports()[..2].to_vec();
    assert!(
        JunctionSchema::try_new(
            StableId::new("bad-skew").unwrap(),
            ports,
            JunctionLaw::Dirac {
                rows: 2,
                incidence: vec![0., 1., 1., 0.]
            }
        )
        .is_err()
    );
}
#[test]
fn rejected_trials_keep_last_accepted_window() {
    let mut a = dynamic();
    a.step_with_force(drive(), 0.001).unwrap();
    let before = a.last_window_transfers().unwrap().to_vec();
    assert!(
        a.step_with_force(
            ApertureDrive {
                upstream_pressure_pa: f64::NAN,
                ..drive()
            },
            0.001
        )
        .is_err()
    );
    assert_eq!(a.last_window_transfers().unwrap(), before);
    let mut c = coupled();
    c.step(&[10., -20., 5.]).unwrap();
    let before = c.last_window_transfers().unwrap().to_vec();
    assert!(c.step(&[f64::NAN, 0., 0.]).is_err());
    assert_eq!(c.last_window_transfers().unwrap(), before);
    let mut p = single();
    p.step(
        piston_fixture::drive(20.),
        piston_fixture::controls(),
        || false,
    )
    .unwrap();
    let before = p.last_window_transfers().unwrap().to_vec();
    assert!(
        p.step(
            piston_fixture::drive(20.),
            piston_fixture::controls(),
            || true
        )
        .is_err()
    );
    assert_eq!(p.last_window_transfers().unwrap(), before);
}
