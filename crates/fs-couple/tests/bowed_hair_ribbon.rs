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
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::MultiBowStringBank;
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

/// Measured horse bow-hair, Pitteroff & Woodhouse 1998 model violin bow.
fn violin_hair() -> HairRibbonProperties {
    HairRibbonProperties {
        youngs_modulus_pa: 7.0e9,
        hair_radius_m: 0.19e-3 / 2.0,
        hairs_in_contact: 50.0,
        free_length_m: 0.65,
        excitation_fraction: 0.5,
        distributed_damping_kg_s_m: 1000.0,
        density_kg_m3: 1100.0,
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
