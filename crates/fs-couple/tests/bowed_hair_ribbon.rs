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
