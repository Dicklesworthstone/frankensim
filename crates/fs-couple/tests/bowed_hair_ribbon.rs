//! Compliant bow-hair ribbon on a material string (IRP pexfp.154, M4.4-W2).
//!
//! A tensioned, compliant hair ribbon (Pitteroff & Woodhouse 1998, Part I) in
//! distributed line contact with the string: each contact cell carries its own
//! friction state driven through a per-cell series spring-dashpot whose
//! stiffness and damping come from the horse bow-hair pack, so the contact
//! velocities differ across the width and differential slip EMERGES. The
//! material/body coefficients here are wiring fixtures; the HAIR numbers are
//! the measured pack values. This file grows the width ladder, point-bow limit
//! and the Pitteroff-Woodhouse trend as those gates land.
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::compliant_ribbon::{
    CompliantBowRibbon, HairRibbonProperties,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringSource, ReciprocalStringBank, StringBankOptions, StringBowLocation,
    plate::{SharedMaterialPlateBody, SharedPlateBridge},
};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::render::schedule::force::coupled::{
    ModalCouplingConfig,
    contact::friction::{FrictionDrive, ModalFrictionConfig, multiple::MultiFrictionConfig},
};
use fs_couple::string_specimen::{
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA, with_uniform_circular_material_state,
};
use fs_couple::thin_plate::{
    PlateMaterialModel, PlateThicknessConstraint, with_uniform_plate_material_state,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    gas::{GasSpec, GasState},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_scenario::{PrestressedString, ThinPlate};
use fs_tribo::FrictionLaw;

const STRING_LENGTH_M: f64 = 0.5;

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 6,
                    kernel_id: 5,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn material(rho: f64) -> ResolvedMaterialStatePoint {
    let (mut claims, mut requirements) = (ClaimSet::new(), Vec::new());
    for (key, dims, value) in [
        ("density", Density::DIMS, rho),
        ("young_modulus", Pressure::DIMS, 200e9),
        ("poisson_ratio", Dims::NONE, 0.3),
        (ETA, DynViscosity::DIMS, 10.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 310.0)
                    .with("omega", 0.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic bow-hair fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn source(rho: f64) -> MaterialStringSource {
    let string = PrestressedString {
        length_m: STRING_LENGTH_M,
        tension_n: 20.0,
        lin_density_kg_m: 1.0,
        axial_stiffness_n: 1.0,
        width_m: 1.0,
        n_modes: 24,
        damping_ratio: 0.0,
        rayleigh: None,
        bending_stiffness_n_m2: 0.0,
        kelvin_voigt_bending: None,
        relaxation_bending: None,
        polarization_detune: 0.0,
        moving_end: false,
        tension_loss_factor: None,
    };
    MaterialStringSource::Homogeneous(
        with_uniform_circular_material_state(string, 0.0003, &material(rho))
            .unwrap()
            .with_kelvin_voigt_bending_loss()
            .unwrap(),
    )
}

/// A velocity-weakening bow law, so the contact genuinely sticks and slips.
fn law() -> FrictionLaw {
    FrictionLaw::Stribeck {
        static_mu: 1.2,
        kinetic_mu: 0.35,
        characteristic_speed: 0.05,
        viscous_per_speed: 0.0,
    }
}

fn options() -> StringBankOptions {
    StringBankOptions {
        sample_rate_hz: 48000,
        bow: StringBowLocation {
            string: 0,
            station_fraction: 0.2,
        },
        string_budget: ModalAcousticTimeBudget::audible_reference(),
        network: ModalCouplingConfig {
            max_modes: 32,
            max_connections: 8,
            max_setup_terms: 16384,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-8,
        },
        friction: ModalFrictionConfig {
            max_iterations: 96,
            maximum_force_n: 100.0,
            maximum_slip_speed_m_s: 10.0,
            force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-9,
            sticking_tolerance_m: 1e-12,
        },
    }
}

fn ribbon_limits() -> MultiFrictionConfig {
    MultiFrictionConfig {
        max_ports: 15,
        max_sweeps: 128,
        max_setup_terms: 1_000_000,
    }
}

/// The Pitteroff & Woodhouse 1998 horse bow-hair matdb pack scalars, keyed by
/// property_key. Same seed-v1 source `bow_hair_pack.rs` gates for physicality;
/// here we READ the measured hair numbers so `violin_hair()` is data, not a
/// hand-typed fixture (pexfp.154 DONE-WHEN 2 provenance; MagentaBadger 6545/6556).
fn bow_hair_pack() -> std::collections::BTreeMap<String, f64> {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("repo root")
        .to_path_buf();
    let text = std::fs::read_to_string(root.join(
        "data/matdb/seed-v1/bow-hair-horse-pitteroff-woodhouse-acustica1998/properties.tsv",
    ))
    .expect("bow-hair properties.tsv");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut m = std::collections::BTreeMap::new();
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        if c[0] == "scalar" {
            m.insert(c[3].to_string(), c[4].parse::<f64>().expect("scalar value"));
        }
    }
    m
}

/// The stated lumped tangential ribbon damping [kg/s] at the stated reference
/// width, for the violin bow (pexfp.154 DONE-WHEN 2, MagentaBadger 6556 item 1).
///
/// P&W explicitly publish D = 10 kg/s on p. 548 and in Table I. On p. 547
/// they derive the single-hair estimate from low-frequency decay tests and
/// state that its frequency-independent continuation is an approximation.
/// This is the paper's reduced model, not a measured broadband viscosity or a
/// matched single-hair dataset. The range-based diagnostic below is only a
/// consistency check; it does not fit or determine the rendered coefficient.
/// Measured horse bow-hair, Pitteroff & Woodhouse 1998 model violin bow, read
/// from the matdb pack. The hair radius is the pack diameter range's MIDPOINT
/// (computed here, MagentaBadger 6556 item 2); gamma = 0.5 is a documented
/// modelling choice (6556 item 3); the damping is the stated assumption above.
fn violin_hair() -> HairRibbonProperties {
    let p = bow_hair_pack();
    // Hair radius = midpoint of the printed diameter range, in code from min/max.
    let diameter_midpoint_m =
        0.5 * (p["hair_diameter_printed_range_low"] + p["hair_diameter_printed_range_high"]);
    let reference_width_m = p["violin_bow_ribbon_width"];
    HairRibbonProperties {
        youngs_modulus_pa: p["young_modulus_static"],
        hair_radius_m: 0.5 * diameter_midpoint_m,
        hairs_in_contact: p["hairs_in_contact_with_string_at_most"],
        free_length_m: p["violin_bow_hair_length_frog_to_tip"],
        excitation_fraction: 0.5, // gamma: a documented modelling choice (6556)
        // Stated assumption D, carried as the distributed d = D / b_ref so the
        // lumped tangential damping D = d * b_ref is width-free (frankensim-o7qbp).
        distributed_damping_kg_s_m: p["nominal_ribbon_damping_pw98_model"] / reference_width_m,
        reference_width_m,
        density_kg_m3: p["density"],
    }
}

/// One material string on a plate body (no cavity), the bow at station 0.2.
fn base(cx: &Cx<'_>, gate: &CancelGate) -> ReciprocalStringBank {
    let plate = with_uniform_plate_material_state(
        ThinPlate {
            length_m: 0.4,
            width_m: 0.3,
            thickness_m: 1.0,
            density_kg_m3: 1.0,
            e1_pa: 1.0,
            e2_pa: 1.0,
            nu12: 0.0,
            g12_pa: 1.0,
            material_angle_rad: 0.0,
            damping_ratio: 0.0,
            thermoelastic: None,
            kelvin_voigt_bending: None,
            n_modes: 2,
            geometric_nonlinearity: false,
            pretension_n_m: 0.0,
            clamped: false,
        },
        &material(7800.0),
        PlateMaterialModel::Isotropic,
        PlateThicknessConstraint::FixedThickness(0.001),
    )
    .unwrap()
    .with_kelvin_voigt_bending_loss()
    .unwrap();
    let body = SharedMaterialPlateBody {
        specimen: plate,
        // One bridge for the one string.
        bridges: vec![SharedPlateBridge {
            station_fraction: 0.8,
            stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1,
            rest_extension_m: 0.0,
        }],
        modal_budget: ModalAcousticTimeBudget::audible_reference(),
        listener_m: 1.0,
    };
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    ReciprocalStringBank::from_material_plate(
        cx,
        vec![source(7800.0)],
        gas,
        body,
        law(),
        options(),
        gate,
    )
    .unwrap()
}

/// Per-cell prescribed drives: common bow speed, total load shared equally.
fn ribbon_drives(cells: usize, speed: f64, load: f64) -> Vec<FrictionDrive> {
    vec![
        FrictionDrive {
            speed_m_s: speed,
            normal_force_n: load / cells as f64
        };
        cells
    ]
}

/// G3: the common performance route must retain the original total bow load,
/// compliant hair state and the released instrument's vibration. Manual joint
/// stepping is an independent check of the schedule-to-physical-port mapping.
#[test]
fn g3_compliant_ribbon_schedules_one_total_load_and_preserves_release() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::MultiBowPerformanceConfig;
    use fs_couple::pcm_wav::observation::PressureRenderer;
    use fs_scenario::gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue};

    with_cx(|cx, gate| {
        let cells = 5;
        let make = || base(cx, gate)
            .with_compliant_bow_ribbon(0.02, cells, violin_hair(), ribbon_limits(), gate)
            .unwrap();
        let mut manual = make();
        let scheduled = make();
        let ribbon = scheduled.source().ribbon().expect("retained contact quadrature");
        assert_eq!(ribbon.stations().len(), cells);
        assert_eq!(ribbon.centre(), options().bow);
        let schedule = GestureSchedule::try_new(1000, vec![GestureTrack {
            id: "stroke".into(),
            target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow {
                velocity_m_per_s: 0.01, normal_force_n: 0.2, station: 0.2,
            },
            events: vec![GestureEvent {
                time_s: 0.004, transition_s: 0.0,
                value: GestureValue::Bow {
                    velocity_m_per_s: 0.0, normal_force_n: 0.0, station: 0.2,
                },
            }],
        }]).unwrap();
        let mut performance = scheduled.into_scheduled_ribbon(schedule, "stroke".into(),
            MultiBowPerformanceConfig { samples: 384, max_block: 48,
                max_work: 100_000, max_events: 1000 }).unwrap();
        let zero = vec![0.0; manual.system().mode_count()];
        let mut expected = Vec::new();
        let mut actual = vec![0.0; 384];
        let mut drive = vec![FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }; cells];
        for sample in 0..384 {
            let total = if sample < 192 {
                FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.2 }
            } else { FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 } };
            manual.source().ribbon().unwrap().project_drive(total, &mut drive).unwrap();
            let frame = manual.step(&zero, &drive).unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            expected.push(frame.observer_pressure_pa);
        }
        for chunk in actual.chunks_mut(48) { performance.block(chunk).unwrap(); }
        assert_eq!(actual, expected);
        assert!(actual[192..].iter().any(|p| p.abs() > 1e-12), "release muted the bank");
        assert_eq!(performance.bank().system().last_frame(), manual.system().last_frame());
    });
}

#[test]
fn g0_compliant_ribbon_runs_with_differential_slip_and_balanced_energy() {
    with_cx(|cx, gate| {
        let (width, cells) = (0.02, 5usize);
        // project_drive on a reference ribbon matches the manual per-cell share.
        let reference =
            CompliantBowRibbon::new(options().bow, STRING_LENGTH_M, width, cells, violin_hair())
                .unwrap();
        let mut projected = vec![
            FrictionDrive {
                speed_m_s: 0.0,
                normal_force_n: 0.0
            };
            cells
        ];
        reference
            .project_drive(
                FrictionDrive {
                    speed_m_s: 0.2,
                    normal_force_n: 1.0,
                },
                &mut projected,
            )
            .unwrap();
        assert_eq!(projected, ribbon_drives(cells, 0.2, 1.0));

        let mut soft = base(cx, gate)
            .with_compliant_bow_ribbon(width, cells, violin_hair(), ribbon_limits(), gate)
            .unwrap();
        // A much stiffer hair, otherwise identical, for the comparison.
        let stiffer = HairRibbonProperties {
            youngs_modulus_pa: 5.0e10,
            ..violin_hair()
        };
        let mut hard = base(cx, gate)
            .with_compliant_bow_ribbon(width, cells, stiffer, ribbon_limits(), gate)
            .unwrap();

        // A gentle, convergent bow on this low-tension wiring string.
        let zero = vec![0.0; soft.system().mode_count()];
        let drives = ribbon_drives(cells, 0.01, 0.2);
        let (mut differential_slip, mut hair_stiffness_matters) = (false, false);
        for _ in 0..600 {
            let s = soft.step(&zero, &drives).unwrap().clone();
            let h = hard.step(&zero, &drives).unwrap().clone();
            // The joint energy ledger still closes with the per-cell hair springs.
            assert!(
                s.energy_residual_j.abs() <= s.energy_tolerance_j,
                "compliant ribbon energy residual {} exceeds tolerance {}",
                s.energy_residual_j,
                s.energy_tolerance_j
            );
            // Differential slip across the width: the cells do not share one slip.
            let (lo, hi) = s
                .ports
                .iter()
                .map(|p| p.slip_distance_m)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                    (lo.min(x), hi.max(x))
                });
            if hi - lo > 1e-12 {
                differential_slip = true;
            }
            // The hair stiffness changes the motion: a soft hair is not a stiff one.
            for i in 0..cells {
                if (s.ports[i].traction_n - h.ports[i].traction_n).abs() > 1e-9 {
                    hair_stiffness_matters = true;
                }
            }
        }
        assert!(
            differential_slip,
            "the finite-width contact showed no differential slip"
        );
        assert!(
            hair_stiffness_matters,
            "hair stiffness did not change the motion"
        );
    });
}

#[test]
fn g0_binding_refusals() {
    with_cx(|cx, gate| {
        // Even cell count and a ribbon wider than the string refuse at construction.
        assert!(
            base(cx, gate)
                .with_compliant_bow_ribbon(0.02, 4, violin_hair(), ribbon_limits(), gate)
                .is_err()
        );
        assert!(
            base(cx, gate)
                .with_compliant_bow_ribbon(0.9, 5, violin_hair(), ribbon_limits(), gate)
                .is_err()
        );
        // Non-physical hair data refuses.
        let mut bad = violin_hair();
        bad.youngs_modulus_pa = 0.0;
        assert!(
            base(cx, gate)
                .with_compliant_bow_ribbon(0.02, 5, bad, ribbon_limits(), gate)
                .is_err()
        );
    });
}

#[test]
fn g0_restore_string_states_seeds_the_motion_and_refuses_bad_input() {
    with_cx(|cx, gate| {
        let (width, cells) = (0.02, 5usize);
        let mut bank = base(cx, gate)
            .with_compliant_bow_ribbon(width, cells, violin_hair(), ribbon_limits(), gate)
            .unwrap();
        let n = bank.system().components()[0].states().len();
        assert!(
            bank.system().total_energy_j().unwrap() <= 1e-12,
            "a fresh bank should start at rest"
        );
        // Seed the string's first mode with a known velocity.
        let mut states = vec![
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: 0.0
            };
            n
        ];
        states[0].velocity_m_sqrt_kg_per_s = 1.0e-3;
        bank.restore_string_states(0, &states).unwrap();
        assert!(
            (bank.system().components()[0].states()[0].velocity_m_sqrt_kg_per_s - 1.0e-3).abs()
                <= 1e-15,
            "the restored velocity did not reach the string's first mode"
        );
        assert!(
            bank.system().total_energy_j().unwrap() > 0.0,
            "restoring a nonzero state must seed motion"
        );
        // A step from the seeded state (released bow, free ring) still balances.
        let zero = vec![0.0; bank.system().mode_count()];
        let released = ribbon_drives(cells, 0.0, 0.0);
        let f = bank.step(&zero, &released).unwrap();
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);

        // Refusals: a wrong mode count and an unknown string before playback.
        let mut fresh = base(cx, gate)
            .with_compliant_bow_ribbon(width, cells, violin_hair(), ribbon_limits(), gate)
            .unwrap();
        assert!(fresh.restore_string_states(0, &states[..n - 1]).is_err());
        assert!(fresh.restore_string_states(9, &states).is_err());
        // And a refusal once a sample has been rendered.
        fresh.step(&zero, &released).unwrap();
        assert!(fresh.restore_string_states(0, &states).is_err());
    });
}

#[test]
fn g0_rank_refusal_when_the_string_underresolves_the_quadrature() {
    // Owner ruling 2026-10-23: a ribbon needs >= 4 string modes per contact
    // cell, else the finite-width joint solve is rank-deficient. The 24-mode
    // fixture string carries 5 cells (needs 20) but refuses 9 (needs 36).
    with_cx(|cx, gate| {
        assert!(
            base(cx, gate)
                .with_compliant_bow_ribbon(0.02, 5, violin_hair(), ribbon_limits(), gate)
                .is_ok()
        );
        assert!(
            base(cx, gate)
                .with_compliant_bow_ribbon(0.02, 9, violin_hair(), ribbon_limits(), gate)
                .is_err()
        );
    });
}

// ---- Width ladder + point-bow limit on a tuned string from a Helmholtz start ----
//
// Owner ruling 2026-10-23: scale the string modes with the quadrature
// (n_modes >= 4 * n_cells per rung) so the finite-width solve is well resolved;
// take the point-bow limit at a carried force; keep the clean Helmholtz duty in
// the .152 fixture. The slip fraction here is a stable carried regime, not
// beta; the gate is its CONVERGENCE with quadrature refinement.

// One low target pitch (owner ruling: gate 1 option 2 at one pitch). At 50 Hz
// the finest reachable rung's top mode (cells 17 -> 68 modes -> 3.4 kHz) stays
// under the material's ~4-7 kHz modal-admission ceiling, so {5,9,17} are all
// admitted at n_modes = 4*cells WITHOUT widening the material (owner ruling).
const PITCH_HZ: f64 = 50.0;
const TUNED_LENGTH_M: f64 = 0.33;
const TUNED_LIN_DENSITY_KG_M: f64 = 5.0e-3;
const BOW_BETA: f64 = 0.11;
const BOW_SPEED_M_S: f64 = 0.1;
const BOW_FORCE_N: f64 = 0.15; // a force the joint solve carries at every rung

fn tuned_tension_n() -> f64 {
    TUNED_LIN_DENSITY_KG_M * (2.0 * TUNED_LENGTH_M * PITCH_HZ).powi(2)
}

fn tuned_source(n_modes: usize) -> MaterialStringSource {
    let string = PrestressedString {
        length_m: TUNED_LENGTH_M,
        tension_n: tuned_tension_n(),
        lin_density_kg_m: TUNED_LIN_DENSITY_KG_M,
        axial_stiffness_n: 1.0,
        width_m: 1.0,
        n_modes,
        damping_ratio: 0.0,
        rayleigh: None,
        bending_stiffness_n_m2: 0.0,
        kelvin_voigt_bending: None,
        relaxation_bending: None,
        polarization_detune: 0.0,
        moving_end: false,
        tension_loss_factor: None,
    };
    MaterialStringSource::Homogeneous(
        with_uniform_circular_material_state(string, 0.0003, &material(7800.0))
            .unwrap()
            .with_kelvin_voigt_bending_loss()
            .unwrap(),
    )
}

fn tuned_options(n_modes: usize) -> StringBankOptions {
    let mut o = options();
    o.bow = StringBowLocation {
        string: 0,
        station_fraction: BOW_BETA,
    };
    // The string modes plus the plate must fit the coupled network budget.
    o.network.max_modes = n_modes + 8;
    o.network.max_setup_terms = 50_000_000;
    o
}

fn tuned_base(
    cx: &Cx<'_>,
    gate: &CancelGate,
    n_modes: usize,
) -> Result<ReciprocalStringBank, String> {
    let plate = with_uniform_plate_material_state(
        ThinPlate {
            length_m: 0.4,
            width_m: 0.3,
            thickness_m: 1.0,
            density_kg_m3: 1.0,
            e1_pa: 1.0,
            e2_pa: 1.0,
            nu12: 0.0,
            g12_pa: 1.0,
            material_angle_rad: 0.0,
            damping_ratio: 0.0,
            thermoelastic: None,
            kelvin_voigt_bending: None,
            n_modes: 2,
            geometric_nonlinearity: false,
            pretension_n_m: 0.0,
            clamped: false,
        },
        &material(7800.0),
        PlateMaterialModel::Isotropic,
        PlateThicknessConstraint::FixedThickness(0.001),
    )
    .unwrap()
    .with_kelvin_voigt_bending_loss()
    .unwrap();
    let body = SharedMaterialPlateBody {
        specimen: plate,
        // A light termination: strong bridge damping degrades the prepared
        // exponential split at the scaled mode counts, and the carried regime
        // only needs to be stable and convergent (the clean Helmholtz duty
        // lives in the .152 fixture per the owner ruling).
        bridges: vec![SharedPlateBridge {
            station_fraction: 0.8,
            stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1,
            rest_extension_m: 0.0,
        }],
        modal_budget: ModalAcousticTimeBudget::audible_reference(),
        listener_m: 1.0,
    };
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    ReciprocalStringBank::from_material_plate(
        cx,
        vec![tuned_source(n_modes)],
        gas,
        body,
        law(),
        tuned_options(n_modes),
        gate,
    )
    .map_err(|e| format!("from_material_plate({n_modes} modes): {e:?}"))
}

/// Ideal Helmholtz modal velocities (the .152 series), q_k = 0.
fn helmholtz_start(n_modes: usize) -> Vec<ModalAcousticState> {
    let scale = (TUNED_LIN_DENSITY_KG_M * TUNED_LENGTH_M / 2.0).sqrt();
    (1..=n_modes)
        .map(|k| {
            let sign = if k % 2 == 1 { 1.0 } else { -1.0 };
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: scale
                    * (2.0 * BOW_SPEED_M_S / (BOW_BETA * core::f64::consts::PI * k as f64))
                    * sign,
            }
        })
        .collect()
}

fn ladder_limits() -> MultiFrictionConfig {
    MultiFrictionConfig {
        max_ports: 40,
        max_sweeps: 128,
        max_setup_terms: 50_000_000,
    }
}

#[test]
fn g1_compute_wall_blocks_the_mode_scaled_ladder() {
    // Gate 1 (DONE-WHEN 1) wants a width ladder {5,9,17,33} cells with the string
    // modes scaled 4x per cell, the finest rung converging to <= 1%. In the
    // material-bank runtime each step is a JOINT many-port friction solve: a
    // bounded sweep (never the 128-sweep cap) over n_modes modes x n_cells ports,
    // run ~10-17 times per step. That cost grows with the product modes x ports,
    // so the full run (three 50 Hz periods, 2880 steps) at the scaled rungs is not
    // computable here. The runtime has TWO walls: the per-step compute cost, and a
    // hard Kelvin-Voigt modal-ADMISSION ceiling that refuses the string at
    // construction once the top scaled mode clears it (and the 33-cell/132-mode
    // rung would also exceed the fs_phs 256 exact-ZOH block cap, bead
    // frankensim-wohg3). This test MEASURES both: per reachable rung it reports
    // sweeps/step, wall-time/step and the projected full-run cost; per refused rung
    // it records the refusal reason. It is the wall measurement the owner ruled is
    // the work here (Option 3 / lightweight runtime refused as gate-1 convergence
    // evidence); the convergence itself is deferred to wohg3. The asserts are
    // machine-independent (a reachable rung never caps; the runtime does run).
    with_cx(|cx, gate| {
        let per_period = (48_000.0 / PITCH_HZ).round() as usize; // 960 at 50 Hz
        let full_run = 3 * per_period; // the ladder's run length, 2880 steps
        // Bound each rung by wall-clock AND step count so the test documents the
        // wall without ever running the (infeasible) full ladder.
        let budget = std::time::Duration::from_secs(6);
        let max_probe_steps = 80usize;
        // (steps, ms/step, mean sweeps, max sweeps, capped steps) or a refusal.
        let probe = |cells: usize| -> Result<(usize, f64, f64, usize, usize), String> {
            let n_modes = 4 * cells;
            let mut bank = tuned_base(cx, gate, n_modes)?
                .with_compliant_bow_ribbon(0.02, cells, violin_hair(), ladder_limits(), gate)
                .map_err(|e| format!("with_compliant_bow_ribbon: {e:?}"))?;
            bank.restore_string_states(0, &helmholtz_start(n_modes))
                .map_err(|e| format!("restore: {e:?}"))?;
            let zero = vec![0.0; bank.system().mode_count()];
            let drives = ribbon_drives(cells, BOW_SPEED_M_S, BOW_FORCE_N);
            let (mut steps, mut sweep_total, mut max_sweeps, mut capped) =
                (0usize, 0usize, 0usize, 0usize);
            let t0 = std::time::Instant::now();
            while steps < max_probe_steps && t0.elapsed() < budget {
                let f = bank.step(&zero, &drives).map_err(|e| format!("step {steps}: {e:?}"))?;
                steps += 1;
                sweep_total += f.sweeps;
                max_sweeps = max_sweeps.max(f.sweeps);
                if f.sweeps >= 128 {
                    capped += 1;
                }
            }
            let ms_per_step = 1.0e3 * t0.elapsed().as_secs_f64() / steps.max(1) as f64;
            let mean_sweeps = sweep_total as f64 / steps.max(1) as f64;
            Ok((steps, ms_per_step, mean_sweeps, max_sweeps, capped))
        };
        let mut reachable = 0usize;
        for &cells in &[5usize, 9, 17] {
            let n_modes = 4 * cells;
            match probe(cells) {
                Ok((steps, ms_per_step, mean_sweeps, max_sweeps, capped)) => {
                    let projected_full_run_s = ms_per_step * full_run as f64 / 1.0e3;
                    println!(
                        "{{\"test\":\"g1_compute_wall\",\"cells\":{cells},\"n_modes\":{n_modes},\"ports\":{cells},\"reached\":true,\"steps\":{steps},\"ms_per_step\":{ms_per_step:.3},\"mean_sweeps\":{mean_sweeps:.2},\"max_sweeps\":{max_sweeps},\"capped\":{capped},\"projected_full_run_s\":{projected_full_run_s:.1}}}"
                    );
                    // A reachable rung's cost is the sweep x modes x ports work, NOT
                    // a capped non-convergence: the joint solve never caps, and it
                    // does make progress (it is slow, not stuck).
                    assert_eq!(capped, 0, "cells {cells}: joint solve hit the 128-sweep cap");
                    assert!(steps > 0, "cells {cells}: no step completed in the budget");
                    reachable += 1;
                }
                Err(reason) => println!(
                    "{{\"test\":\"g1_compute_wall\",\"cells\":{cells},\"n_modes\":{n_modes},\"ports\":{cells},\"reached\":false,\"reason\":{reason:?}}}"
                ),
            }
        }
        // The material bank DOES run the finite-width ribbon (so the mechanism is
        // sound), but the mode-scaled ladder hits the two walls above; at least the
        // coarsest rung must be reachable for the measurement to mean anything.
        assert!(reachable >= 1, "no rung was reachable: the runtime did not run at all");
    });
}

/// o7qbp fix: the 5-cell narrow compliant ribbon must now step through the first
/// stick<->slip transition (step 106) instead of refusing, because the active-set
/// rescue resolves the joint partition the Gauss-Seidel sweep cannot converge.
/// RED before the fix (refuses at 106); GREEN after (converges, and the rescue
/// fires at least once).
#[test]
fn narrow_ribbon_converges_through_the_first_transition_via_active_set() {
    with_cx(|cx, gate| {
        let mut bank = tuned_base(cx, gate, 20)
            .unwrap()
            .with_compliant_bow_ribbon(0.0008, 5, violin_hair(), ladder_limits(), gate)
            .unwrap();
        bank.restore_string_states(0, &helmholtz_start(20)).unwrap();
        let zero = vec![0.0; bank.system().mode_count()];
        let drives = ribbon_drives(5, BOW_SPEED_M_S, BOW_FORCE_N);
        let mut total_rescues = 0usize;
        for s in 0..300 {
            let f = bank
                .step(&zero, &drives)
                .unwrap_or_else(|e| panic!("narrow refused at step {s}: {e:?}"));
            total_rescues += f.active_set_rescues;
        }
        println!("{{\"test\":\"narrow_active_set\",\"steps\":300,\"total_rescues\":{total_rescues}}}");
        assert!(
            total_rescues > 0,
            "active-set rescue never fired; the GS sweep converged alone (did the o7qbp repro regress?)"
        );
    });
}

/// o7qbp item 1-2: when the active-set pivot budget is too small to reach the
/// unique consistent partition, the step refuses with the NAMED active-set error
/// (not the original Gauss-Seidel error). The default budget converges.
#[test]
fn active_set_refuses_with_a_named_error_when_the_budget_is_too_small() {
    with_cx(|cx, gate| {
        let drives = ribbon_drives(5, BOW_SPEED_M_S, BOW_FORCE_N);
        // Default budget: the narrow converges through the first transition.
        let mut ok = tuned_base(cx, gate, 20)
            .unwrap()
            .with_compliant_bow_ribbon(0.0008, 5, violin_hair(), ladder_limits(), gate)
            .unwrap();
        ok.restore_string_states(0, &helmholtz_start(20)).unwrap();
        let zero = vec![0.0; ok.system().mode_count()];
        for s in 0..150 {
            ok.step(&zero, &drives)
                .unwrap_or_else(|e| panic!("default budget refused at step {s}: {e:?}"));
        }
        // Budget 1: the rescue cannot resolve the partition, so the step refuses
        // with the named active-set error.
        let mut starved = tuned_base(cx, gate, 20)
            .unwrap()
            .with_compliant_bow_ribbon(0.0008, 5, violin_hair(), ladder_limits(), gate)
            .unwrap();
        starved.force_active_set_budget(1);
        starved.restore_string_states(0, &helmholtz_start(20)).unwrap();
        let mut refusal = None;
        for s in 0..150 {
            if let Err(e) = starved.step(&zero, &drives) {
                refusal = Some((s, format!("{e:?}")));
                break;
            }
        }
        let (s, msg) = refusal.expect("budget-1 narrow must refuse at the first transition");
        println!("{{\"test\":\"budget1\",\"refused_at\":{s},\"msg\":{msg:?}}}");
        // The refusal must be the NAMED active-set exhaustion carrying its payload
        // (pivots against budget, last active set, worst KKT violation), never the
        // original Gauss-Seidel error.
        assert!(
            msg.contains("ActiveSetExhausted")
                && msg.contains("pivots")
                && msg.contains("budget")
                && msg.contains("last_active_set")
                && msg.contains("worst_port")
                && msg.contains("worst_condition"),
            "refusal must be the named ActiveSetExhausted with its payload, got: {msg}"
        );
    });
}

/// o7qbp item 3 (bit-identity): a step the Gauss-Seidel sweep converges never
/// invokes the active-set rescue, so the fix cannot perturb it. Run the
/// converging 1-cell point bow twice -- once with the rescue at its default
/// budget, once with the budget forced to 0 (rescue disabled) -- and require the
/// two traction traces byte-for-byte identical. The rescue-count is 0 on every
/// step, proving the active-set path is inert on a converged run; the trace hash
/// is reported for the record.
#[test]
fn active_set_is_inert_on_a_gs_converged_run() {
    with_cx(|cx, gate| {
        let trace = |budget: Option<usize>| -> Vec<u64> {
            let mut bank = tuned_base(cx, gate, 20)
                .unwrap()
                .with_compliant_bow_ribbon(0.0008, 1, violin_hair(), ladder_limits(), gate)
                .unwrap();
            if let Some(b) = budget {
                bank.force_active_set_budget(b);
            }
            bank.restore_string_states(0, &helmholtz_start(20)).unwrap();
            let zero = vec![0.0; bank.system().mode_count()];
            let drives = ribbon_drives(1, BOW_SPEED_M_S, BOW_FORCE_N);
            let mut bits = Vec::with_capacity(200);
            for s in 0..200 {
                let f = bank
                    .step(&zero, &drives)
                    .unwrap_or_else(|e| panic!("point bow refused at step {s}: {e:?}"));
                assert_eq!(
                    f.active_set_rescues, 0,
                    "a converged GS step must not invoke the active-set rescue (step {s})"
                );
                bits.push(f.ports[0].traction_n.to_bits());
            }
            bits
        };
        let with_rescue = trace(None); // default budget: rescue available
        let rescue_off = trace(Some(0)); // budget 0: rescue disabled
        assert_eq!(
            with_rescue, rescue_off,
            "the active-set changed a Gauss-Seidel-converged run; it must be inert there"
        );
        let hash = with_rescue
            .iter()
            .enumerate()
            .fold(0u64, |h, (s, &b)| h ^ b.rotate_left((s % 64) as u32));
        println!("{{\"test\":\"bit_identity\",\"point_bow_traction_hash\":\"{hash:#018x}\"}}");
    });
}

// ---- Gate 2 (DONE-WHEN 2): point-bow limit on SMOOTH observables ----
//
// MagentaBadger 6545/6556. The quantized port-0 sliding fraction is a FINDING,
// not the gate: its 1/tail step gives false non-monotone "orders". The gate is
// the width -> 0 point-bow limit on SMOOTH observables, against the width-free
// 1-cell compliant point reference:
//   (i)  tail-mean TOTAL bow traction [N] = sum of the port tractions, the net
//        tangential force on the string. Port 0 alone is NOT cell-count
//        invariant: ribbon_drives shares the normal load as load/cells, so a
//        5-cell port 0 carries 1/5 of a 1-cell port 0. For coincident ports
//        (w -> 0) the cells each slide at mu(v)*N/cells and sum to mu(v)*N, the
//        1-cell point's single-port traction, so the TOTAL converges to the
//        reference and equals port 0 for the 1-cell reference itself.
//   (ii) the steady-state period of the observer pressure (the body/bridge
//        signal the frame exposes), by autocorrelation refined to sub-sample.
// PASS = both within 2% of the reference on the TWO FINEST rungs, compared only
// within the same motion regime. The measured order is reported, not asserted;
// monotonicity is not asserted. Sliding fraction (with its 1/tail quantization)
// and the regime are reported per rung.

const SAMPLE_RATE_HZ: f64 = 48_000.0;
const SLIP_EPS_M: f64 = 1.0e-12; // |per-step slip| above this is sliding (travel/step ~2e-6 m)
// One FIXED string resolution for the reference and every rung, so the only
// quantity that changes down the ladder is the contact width. 20 = 4 * 5 cells
// is the owner's minimum (n_modes >= 4*n_cells) for the finest 5-cell ribbon.
const GATE_MODES: usize = 20;

#[derive(Clone, Debug)]
struct RungObs {
    width_m: f64,
    tail_mean_total_traction_n: f64,
    tail_mean_port0_traction_n: f64,
    period_s: f64,
    sliding_fraction_port0: f64,
    slips_per_period_port0: f64,
    regime: &'static str,
    rescues: usize,
}

/// Dominant period [s] of a tail signal by autocorrelation around the expected
/// fundamental, refined to sub-sample by a parabolic fit at the peak. The search
/// window 0.5x..2.5x the expected period admits period-doubling (multiple slip).
fn dominant_period_s(signal: &[f64], expected_period_samples: f64) -> f64 {
    let n = signal.len();
    if n < 4 {
        return 0.0;
    }
    let mean = signal.iter().sum::<f64>() / n as f64;
    let x: Vec<f64> = signal.iter().map(|v| v - mean).collect();
    let acf = |lag: usize| -> f64 {
        let mut s = 0.0;
        for i in 0..(n - lag) {
            s += x[i] * x[i + lag];
        }
        s
    };
    let lo = ((0.5 * expected_period_samples).round() as usize).max(1);
    let hi = ((2.5 * expected_period_samples).round() as usize).min(n - 2);
    let mut best_lag = lo;
    let mut best = f64::NEG_INFINITY;
    for lag in lo..=hi {
        let r = acf(lag);
        if r > best {
            best = r;
            best_lag = lag;
        }
    }
    // Parabolic sub-sample refinement at the peak.
    let (rm, r0, rp) = (acf(best_lag - 1), acf(best_lag), acf(best_lag + 1));
    let denom = rm - 2.0 * r0 + rp;
    let delta = if denom.abs() > 0.0 {
        0.5 * (rm - rp) / denom
    } else {
        0.0
    };
    (best_lag as f64 + delta) / SAMPLE_RATE_HZ
}

fn classify_regime(slips_per_period: f64, period_s: f64, expected_period_samples: f64) -> &'static str {
    let period_ratio = period_s * SAMPLE_RATE_HZ / expected_period_samples; // 1.0 = fundamental
    if slips_per_period < 0.5 {
        "stuck_or_aperiodic"
    } else if (0.7..=1.4).contains(&slips_per_period) && (0.8..=1.25).contains(&period_ratio) {
        "helmholtz_one_slip"
    } else if slips_per_period <= 3.5 {
        "multiple_slip"
    } else {
        "raucous_or_aperiodic"
    }
}

/// Run one Helmholtz-started rung and measure the smooth observables over the
/// tail. `n_modes` is FIXED across the reference and every rung so the only
/// quantity that varies is the contact width (not the string resolution); it
/// must satisfy the owner rule n_modes >= 4*cells for the finest ribbon.
fn run_rung(
    cx: &Cx<'_>,
    gate: &CancelGate,
    hair: HairRibbonProperties,
    width_m: f64,
    cells: usize,
    n_modes: usize,
    run: usize,
    tail: usize,
) -> Result<RungObs, String> {
    let mut bank = tuned_base(cx, gate, n_modes)?
        .with_compliant_bow_ribbon(width_m, cells, hair, ladder_limits(), gate)
        .map_err(|e| format!("with_compliant_bow_ribbon(w={width_m:e}, cells={cells}): {e:?}"))?;
    bank.restore_string_states(0, &helmholtz_start(n_modes))
        .map_err(|e| format!("restore: {e:?}"))?;
    let zero = vec![0.0; bank.system().mode_count()];
    let drives = ribbon_drives(cells, BOW_SPEED_M_S, BOW_FORCE_N);
    let tail_start = run - tail;
    let (mut total_tr, mut port0_tr, mut pressure) = (
        Vec::with_capacity(tail),
        Vec::with_capacity(tail),
        Vec::with_capacity(tail),
    );
    let (mut sliding, mut slip_events, mut prev_sliding, mut rescues) = (0usize, 0usize, false, 0usize);
    for s in 0..run {
        let f = bank.step(&zero, &drives).map_err(|e| format!("step {s}: {e:?}"))?;
        rescues += f.active_set_rescues;
        if s >= tail_start {
            total_tr.push(f.ports.iter().map(|p| p.traction_n).sum::<f64>());
            port0_tr.push(f.ports[0].traction_n);
            pressure.push(f.observer_pressure_pa);
            let is_sliding = f.ports[0].slip_distance_m.abs() > SLIP_EPS_M;
            if is_sliding {
                sliding += 1;
            }
            if is_sliding && !prev_sliding {
                slip_events += 1;
            }
            prev_sliding = is_sliding;
        }
    }
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len().max(1) as f64;
    let expected_period = SAMPLE_RATE_HZ / PITCH_HZ;
    let period_s = dominant_period_s(&pressure, expected_period);
    let n_periods = tail as f64 / expected_period;
    let slips_per_period = slip_events as f64 / n_periods;
    let sliding_fraction = sliding as f64 / tail as f64;
    Ok(RungObs {
        width_m,
        tail_mean_total_traction_n: mean(&total_tr),
        tail_mean_port0_traction_n: mean(&port0_tr),
        period_s,
        sliding_fraction_port0: sliding_fraction,
        slips_per_period_port0: slips_per_period,
        regime: classify_regime(slips_per_period, period_s, expected_period),
        rescues,
    })
}

/// Print the damping provenance cross-check (MagentaBadger 6556 item 1): the
/// stated value, the eigenfrequency-consistent 2 M alpha estimate, the full
/// (M, alpha) bracket, and S at the diameter range ends (6556 item 2). All from
/// the pack; no gate run. The gate uses the published model's retained source row.
fn print_damping_provenance(hair: &HairRibbonProperties) {
    let p = bow_hair_pack();
    let n_h = p["hairs_in_contact_with_string_at_most"];
    let (al_lo, al_hi) = (
        p["exponential_decay_rate_printed_range_low"],
        p["exponential_decay_rate_printed_range_high"],
    );
    // M and f0 ranges are stated in the pw98-dynamic observation (not scalar
    // rows): attached mass 11..62 g, eigenfrequency 10..100 Hz.
    let (m_lo, m_hi) = (0.011_f64, 0.062_f64);
    let (f_lo, f_hi) = (10.0_f64, 100.0_f64);
    let gamma = hair.excitation_fraction;
    // Lumped tangential stiffness S = n_H pi a^2 E / (l gamma (1-gamma)) [N/m]
    // (same formula as the private HairRibbonProperties::lumped_tangential_stiffness_n_m).
    let s_at = |a: f64| {
        n_h * core::f64::consts::PI * a * a * hair.youngs_modulus_pa
            / (hair.free_length_m * gamma * (1.0 - gamma))
    };
    let s = s_at(hair.hair_radius_m);
    let (al, m, f0) = (0.5 * (al_lo + al_hi), 0.5 * (m_lo + m_hi), 0.5 * (f_lo + f_hi));
    let w0 = 2.0 * core::f64::consts::PI * f0;
    let d_eigen = 2.0 * al * s / (w0 * w0); // M from eigenfrequency, = Kelvin-Voigt eta S / w
    let d_mid_product = n_h * 2.0 * m * al; // unmatched midpoint product
    let bracket = [
        n_h * 2.0 * m_lo * al_lo,
        n_h * 2.0 * m_lo * al_hi,
        n_h * 2.0 * m_hi * al_lo,
        n_h * 2.0 * m_hi * al_hi,
    ];
    let a_lo = 0.5 * p["hair_diameter_printed_range_low"];
    let a_hi = 0.5 * p["hair_diameter_printed_range_high"];
    println!(
        "{{\"provenance\":\"ribbon_damping\",\"stated_D_kg_s\":{},\"derived_D_eigen_consistent_kg_s\":{:.3},\"midpoint_product_D_kg_s\":{:.3},\"bracket_kg_s\":[{:.3},{:.3}],\"S_n_m\":{:.1},\"S_at_diameter_lo_n_m\":{:.1},\"S_at_diameter_hi_n_m\":{:.1},\"hair_radius_m\":{:e}}}",
        p["nominal_ribbon_damping_pw98_model"],
        d_eigen,
        d_mid_product,
        bracket.iter().cloned().fold(f64::INFINITY, f64::min),
        bracket.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        s,
        s_at(a_lo),
        s_at(a_hi),
        hair.hair_radius_m,
    );
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(f64::MIN_POSITIVE)
}

/// Gate 2: the finite-width ribbon's SMOOTH observables converge to the
/// width-free 1-cell compliant point as the width goes to zero; the limit holds
/// on the two finest rungs (within 2% on total traction AND period).
#[test]
fn g2_point_bow_limit_on_smooth_observables() {
    with_cx(|cx, gate| {
        let per_period = (SAMPLE_RATE_HZ / PITCH_HZ).round() as usize; // 960
        let (run, tail) = (4 * per_period, 3 * per_period); // settle 1, measure 3 periods
        print_damping_provenance(&violin_hair());
        // Width-free reference: the 1-cell compliant point bow.
        let reference = run_rung(cx, gate, violin_hair(), 8.0e-4, 1, GATE_MODES, run, tail)
            .unwrap_or_else(|e| panic!("reference point bow refused: {e}"));
        println!(
            "{{\"test\":\"g2_reference_point\",\"total_traction_n\":{},\"port0_traction_n\":{},\"period_s\":{},\"sliding_fraction\":{},\"slips_per_period\":{:.3},\"regime\":{:?}}}",
            reference.tail_mean_total_traction_n,
            reference.tail_mean_port0_traction_n,
            reference.period_s,
            reference.sliding_fraction_port0,
            reference.slips_per_period_port0,
            reference.regime,
        );
        // Width ladder (5 cells) toward the point-bow limit; two finest 7e-5, 5e-5.
        let widths = [8.0e-4_f64, 4.0e-4, 2.0e-4, 1.0e-4, 7.0e-5, 5.0e-5];
        let quantization_floor = 1.0 / tail as f64;
        let mut rungs: Vec<RungObs> = Vec::new();
        for &w in &widths {
            let r = run_rung(cx, gate, violin_hair(), w, 5, GATE_MODES, run, tail)
                .unwrap_or_else(|e| panic!("rung w={w:e} refused: {e}"));
            println!(
                "{{\"test\":\"g2_rung\",\"w\":{:e},\"cells\":5,\"total_traction_n\":{},\"rel_traction\":{:.6},\"port0_traction_n\":{},\"period_s\":{},\"rel_period\":{:.6},\"sliding_fraction\":{},\"quantization_floor\":{:.6},\"slips_per_period\":{:.3},\"regime\":{:?},\"rescues\":{}}}",
                w,
                r.tail_mean_total_traction_n,
                rel(r.tail_mean_total_traction_n, reference.tail_mean_total_traction_n),
                r.tail_mean_port0_traction_n,
                r.period_s,
                rel(r.period_s, reference.period_s),
                r.sliding_fraction_port0,
                quantization_floor,
                r.slips_per_period_port0,
                r.regime,
                r.rescues,
            );
            rungs.push(r);
        }
        // Measured order per halving on the total traction gap (reported only).
        let orders: Vec<f64> = (1..rungs.len())
            .map(|i| {
                let g0 = rel(rungs[i - 1].tail_mean_total_traction_n, reference.tail_mean_total_traction_n);
                let g1 = rel(rungs[i].tail_mean_total_traction_n, reference.tail_mean_total_traction_n);
                let ratio = widths[i - 1] / widths[i];
                if g1 > 0.0 && ratio > 1.0 { (g0 / g1).ln() / ratio.ln() } else { f64::NAN }
            })
            .collect();
        println!("{{\"test\":\"g2_measured_order_traction\",\"orders\":{orders:?}}}");
        // PASS: the two finest rungs within 2% of the reference on BOTH smooth
        // observables, compared within the same regime.
        for r in &rungs[rungs.len() - 2..] {
            assert_eq!(
                r.regime, reference.regime,
                "rung w={:e} regime {:?} differs from reference {:?}; compare only within a regime",
                r.width_m, r.regime, reference.regime
            );
            let rt = rel(r.tail_mean_total_traction_n, reference.tail_mean_total_traction_n);
            let rp = rel(r.period_s, reference.period_s);
            assert!(
                rt <= 0.02,
                "rung w={:e}: total traction {} off reference {} by {:.3}% > 2%",
                r.width_m, r.tail_mean_total_traction_n, reference.tail_mean_total_traction_n, rt * 100.0
            );
            assert!(
                rp <= 0.02,
                "rung w={:e}: period {} off reference {} by {:.3}% > 2%",
                r.width_m, r.period_s, reference.period_s, rp * 100.0
            );
        }
    });
}
