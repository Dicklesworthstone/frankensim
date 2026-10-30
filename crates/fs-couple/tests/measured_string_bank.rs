//! Measured whole strings use the original gas, reciprocal-body and work owners.
//! Both the dry section and receiver are explicitly synthetic test fixtures.
//! This does not de-embed published in-air string fits or model a cello body.
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::bowed_string::BowedStringCard;
use fs_couple::bowed_string::runtime::reciprocal::InteriorStringConnection;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections;
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringMember, MaterialStringSource, ReciprocalStringBank, StringBankOptions,
    StringBowLocation,
};
use fs_couple::material_assembly::measured_string::{
    MeasuredStringProperties, ResolvedMeasuredStringSpecimen,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, ModalFrictionConfig,
};
use fs_couple::render::schedule::force::coupled::{
    ModalAttachment, ModalCouplingConfig, PreparedModalBudget,
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};
use fs_math::c64::C64;
use fs_tribo::FrictionLaw;

fn scoped<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 10,
                    kernel_id: 306,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}
fn measured(n_modes: usize) -> ResolvedMeasuredStringSpecimen {
    let properties = MeasuredStringProperties {
        length: "L".into(),
        diameter: "d".into(),
        tension: "T".into(),
        linear_density: "mu".into(),
        bending_stiffness: "EI".into(),
    };
    let requirements = properties.requirements().unwrap();
    let mut claims = ClaimSet::new();
    // These numbers coincide with the Dominant table but define an authored
    // dry section here. The actual Lampis pack has effective fitted mu/EI and
    // an unknown atmosphere; it is not qualified as a dry fluid-loading input.
    for (requirement, value) in requirements
        .iter()
        .zip([0.70, 0.00115, 116.3, 0.0062, 0.000303])
    {
        let dims = requirement.dims();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(requirement.name(), dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("temperature_known", 0.0, 0.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source: "synthetic dry-section integration fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
                observations: vec![],
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dry uniform section".into(),
            phase: "solid composite".into(),
            process: "authored dry coefficients".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let state = resolve_material_state_point(
        &card,
        &QueryPoint::new().with("temperature_known", 0.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedMeasuredStringSpecimen::try_new(n_modes, &state, properties).unwrap()
}
fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}
fn build(cx: &Cx<'_>, gate: &CancelGate, medium: GasState) -> ReciprocalStringBank {
    let body = ModalAcousticTimeModel::try_new(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: core::f64::consts::TAU * 380.0,
            damping_ratio: 0.01,
            pressure_per_modal_velocity: C64::new(1.0, 0.0),
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let strings = [0.78, 0.83].map(|station_fraction| MaterialStringMember {
        source: MaterialStringSource::Measured(measured(1))
            .with_transverse_polarizations(
                TransverseDirections::try_new([1.0, 0.0], [0.6, 0.8]).unwrap(),
            )
            .unwrap(),
        connections: vec![InteriorStringConnection {
            station_fraction,
            body: ModalAttachment {
                component: 0,
                shapes: vec![1.0 / 0.2_f64.sqrt()],
            },
            stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1,
            rest_extension_m: 0.0,
        }],
    });
    let options = StringBankOptions {
        sample_rate_hz: 48_000,
        bow: StringBowLocation {
            string: 0,
            station_fraction: 0.2,
        },
        string_budget: ModalAcousticTimeBudget::audible_reference(),
        network: ModalCouplingConfig {
            max_modes: 64,
            max_connections: 64,
            max_setup_terms: 1_000_000,
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
    };
    ReciprocalStringBank::new(
        cx,
        strings.to_vec(),
        medium,
        vec![body],
        FrictionLaw::Coulomb {
            static_mu: 0.8,
            kinetic_mu: 0.4,
        },
        options,
        gate,
    )
    .unwrap()
    .with_cylinder_gas_loading(
        CylinderFlowResolution {
            minimum_wavenumber_radius: 1e-3,
            maximum_wavenumber_radius: 30.0,
            intervals: 2,
        },
        gate,
    )
    .unwrap()
    .with_prepared_linear_propagation(
        PreparedModalBudget {
            max_modes: 64,
            max_dense_bytes: 67_108_864,
            max_multiply_terms: 200_000_000_000,
        },
        gate,
    )
    .unwrap()
}

#[test]
fn measured_dispersion_and_actual_radius_reach_sampled_and_causal_gas_paths() {
    scoped(|cx, gate| {
        let source = measured(2);
        let dry = source.string();
        let sampled =
            BowedStringCard::from_measured_string_specimen(cx, &source, &gas(101_325.0), 48_000)
                .unwrap();
        assert!(sampled.zetas.iter().all(|zeta| *zeta > 0.0));
        assert_eq!(sampled.bending_stiffness_n_m2, dry.bending_stiffness_n_m2);
        assert_eq!(sampled.viscous_bending_n_m2_s, 0.0);
        let a = build(cx, gate, gas(101_325.0));
        let b = build(cx, gate, gas(4.0 * 101_325.0));
        let added = core::f64::consts::PI * source.radius_m().powi(2) * gas(101_325.0).density;
        for (bank, factor) in [(&a, 1.0), (&b, 4.0)] {
            let loading = &bank.source().cylinder_gas_loading().unwrap().strings()[0];
            assert_eq!(
                loading.solid_card().linear_density_kg_m,
                dry.lin_density_kg_m
            );
            assert_eq!(loading.flow().radius_m(), source.radius_m());
            let card = &bank.cards()[0];
            assert!(
                (card.linear_density_kg_m / (dry.lin_density_kg_m + factor * added) - 1.0).abs()
                    < 2e-15
            );
            assert_eq!(card.tension_n, dry.tension_n);
            assert_eq!(card.bending_stiffness_n_m2, dry.bending_stiffness_n_m2);
            assert_eq!(
                bank.system().components()[0].modes()[0],
                bank.system().components()[0].modes()[1]
            );
        }
        let ratio = ((dry.lin_density_kg_m + added) / (dry.lin_density_kg_m + 4.0 * added)).sqrt();
        assert!(
            (b.cards()[0].fundamental_hz() / a.cards()[0].fundamental_hz() - ratio).abs() < 2e-14
        );
    });
}

#[test]
fn measured_string_force_reaches_sympathetic_string_fluid_and_body_with_closed_work() {
    scoped(|cx, gate| {
        let mut a = build(cx, gate, gas(101_325.0));
        let mut b = build(cx, gate, gas(4.0 * 101_325.0));
        let mut force = vec![0.0; a.system().mode_count()];
        let released = FrictionDrive {
            speed_m_s: 0.0,
            normal_force_n: 0.0,
        };
        let mut work = 0.0;
        let mut dissipated = 0.0;
        let mut sympathetic = 0.0_f64;
        let mut fluid = 0.0_f64;
        let mut pressure_change = 0.0_f64;
        for step in 0..32 {
            force[0] = if step < 8 { 1.0 } else { 0.0 };
            let fa = *a.step(&force, released).unwrap();
            let fb = *b.step(&force, released).unwrap();
            work += fa.external_work_j;
            dissipated += fa.network_dissipation_j + fa.friction_dissipation_j;
            assert_eq!(fa.drive_work_j, 0.0);
            assert!(fa.energy_residual_j.abs() <= fa.energy_tolerance_j);
            assert!(fb.energy_residual_j.abs() <= fb.energy_tolerance_j);
            sympathetic = sympathetic.max(
                a.system().components()[1]
                    .states()
                    .iter()
                    .map(|s| s.velocity_m_sqrt_kg_per_s.abs())
                    .fold(0.0, f64::max),
            );
            fluid = fluid.max(
                a.system().components()[3..]
                    .iter()
                    .flat_map(|c| c.states())
                    .map(|s| 0.5 * s.velocity_m_sqrt_kg_per_s.powi(2))
                    .sum::<f64>(),
            );
            pressure_change =
                pressure_change.max((fa.observer_pressure_pa - fb.observer_pressure_pa).abs());
        }
        assert!(sympathetic > 0.0 && fluid > 0.0 && pressure_change > 1e-18);
        assert!(work > 1e-9 && dissipated > 0.0);
        assert!(
            (a.system().total_energy_j().unwrap() + dissipated - work).abs() < 2e-13 + 1e-8 * work
        );
    });
}

#[test]
fn measured_transverse_source_refuses_unmeasured_shear_and_axial_mechanics() {
    let source = MaterialStringSource::Measured(measured(2));
    assert!(source.clone().with_torsion("G", "eta_s", 1).is_err());
    scoped(|cx, gate| {
        assert!(
            source
                .into_nonlinear(
                    cx,
                    gas(101_325.0),
                    48_000,
                    0.1,
                    vec![
                        ModalAcousticState {
                            displacement_m_sqrt_kg: 0.0,
                            velocity_m_sqrt_kg_per_s: 0.0
                        };
                        2
                    ],
                    gate
                )
                .is_err()
        );
        gate.request();
        assert!(
            BowedStringCard::from_measured_string_specimen(
                cx,
                &measured(2),
                &gas(101_325.0),
                48_000
            )
            .is_err()
        );
    });
}
