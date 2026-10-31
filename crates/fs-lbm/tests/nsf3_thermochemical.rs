//! Direct behavior regressions of the existing NSF solver with the original
//! named-mixture law and manufactured NASA9/transport source data. These are
//! integration and conservation checks, not experimental gas validation.
#[path = "nsf3_thermochemical/source.rs"]
mod source;

use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{
    NamedMixtureNsfGas, NsfError, NsfGasMedium, NsfGrid3, NsfOpenState, NsfSgsConfig, NsfSpec,
    ThermochemicalNsfGas, WallThermal,
};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse};
use fs_material::gas::GasSpec;
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_qty::semantic::{Composition, CompositionBasis};

fn controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 600.0,
        energy_tolerance_j_kg: 1e-6,
        maximum_iterations: 64,
    }
}
fn gas(offset: f64, pressure_max: f64, k_scale: f64) -> ThermochemicalNsfGas {
    ThermochemicalNsfGas::new(
        source::model(offset, pressure_max, k_scale),
        source::point(300.0, 100000.0),
        controls(),
    )
    .unwrap()
}
fn spec<G>(gas: G) -> NsfSpec<G> {
    NsfSpec {
        dims: [8, 2, 2],
        periodic: [true; 3],
        gas,
        temperature_k: 300.0,
        pressure_pa: 100000.0,
        dx_m: 5e-5,
        temperature_window_k: [250.0, 600.0],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 1.0,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored synthetic ratio; not empirical",
        maximum_cells: 32,
        wall: WallThermal::Adiabatic,
        workers: 1,
        transport_points: 5,
    }
}
fn varying(x: [usize; 3]) -> (f64, [f64; 3], f64) {
    let phase = 2.0 * core::f64::consts::PI * x[0] as f64 / 8.0;
    (
        1.0 + 0.01 * phase.cos(),
        [0.02, 0.0, 0.0],
        300.0 + 0.2 * phase.sin(),
    )
}

#[test]
fn g0_constructor_refuses_a_caloric_bracket_beyond_source_support() {
    let model = source::model(-4000.0, 250000.0, 1.0);
    let point = source::point(300.0, 100000.0);
    let admitted = GasEnergySolve { maximum_temperature_k: 700.0, ..controls() };
    // The reference is supported in both cases, but only the first complete
    // run window is supported. Refuse the other before constructing a grid.
    assert!(ThermochemicalNsfGas::new(model.clone(), point.clone(), admitted).is_ok());
    let unsupported = GasEnergySolve {
        maximum_temperature_k: 700.0_f64.next_up(),
        ..admitted
    };
    assert!(matches!(
        ThermochemicalNsfGas::new(model, point, unsupported),
        Err(NsfError::Gas(_))
    ));
}

#[test]
fn g0_open_faces_use_source_caloric_energy_and_pressure_support() {
    // Uniform flow at a temperature away from the reference detects a ghost
    // energy reconstructed with frozen cv instead of the source's integral.
    // Arbitrary formation energy must not leak into the open energy flux.
    for offset in [-4000.0, 1_000_000.0] {
        let mut s = spec(gas(offset, 110000.0, 1.0));
        s.periodic[0] = false;
        let mut grid =
            NsfGrid3::<D3q27>::new(s, |_| false, |_| (300.0 / 310.0, [0.02, 0.0, 0.0], 310.0))
                .unwrap();
        let exterior = NsfOpenState {
            density_kg_m3: grid.law().rho_ref_kg_m3 * 300.0 / 310.0,
            velocity_m_s: [0.02 * grid.law().dx_m / grid.law().dt_s, 0.0, 0.0],
            temperature_k: 310.0,
        };
        for high in [false, true] {
            grid.set_open_face(0, high, Some(exterior)).unwrap();
        }
        let hash = grid.state_hash();
        let unsupported = NsfOpenState {
            density_kg_m3: 2.0 * exterior.density_kg_m3,
            ..exterior
        };
        assert!(matches!(
            grid.set_open_face(0, false, Some(unsupported)),
            Err(NsfError::Gas(_))
        ));
        assert_eq!(grid.state_hash(), hash);
        for _ in 0..32 {
            let ledger = grid.step().unwrap();
            assert!(ledger.residual_rel < 1e-12);
        }
        for x in 0..8 {
            let (rho, u, temperature) = grid.cell([x, 0, 0]);
            assert!((rho - 300.0 / 310.0).abs() < 1e-10);
            assert!((u[0] - 0.02).abs() < 1e-10);
            assert!((temperature - 310.0).abs() < 1e-6);
        }
    }
}

#[test]
fn nasa9_energy_heating_updates_temperature_pressure_gamma_and_transport_together() {
    let gas = gas(-4000.0, 250000.0, 1.0);
    let reference = gas.model().evaluate(gas.context()).unwrap();
    assert!(reference.specific_internal_energy_j_kg() < 0.0);
    let r = reference.state().specific_gas_constant;
    let dx: f64 = 5e-5;
    let dt = dx / (3.0 * r * 300.0).sqrt();
    let specific_unit = (dx / dt).powi(2);
    let heat_per_mass = 10000.0;
    let source_e = heat_per_mass / specific_unit;
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        spec(gas),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(move |_, _| (0.0, [0.0; 3], source_e))),
    )
    .unwrap();
    let before = grid.total_energy_j();
    let mass = 32.0 * grid.law().rho_ref_kg_m3 * dx.powi(3);
    let steps = 4;
    for _ in 0..steps {
        grid.step().unwrap();
    }
    let (rho, _, temperature) = grid.cell([0, 0, 0]);
    // Independent integral of cv/R = 2 + .001 T, not another runtime EOS.
    let change =
        r * (2.0 * (temperature - 300.0) + 0.0005 * (temperature * temperature - 300.0 * 300.0));
    assert!((change - steps as f64 * heat_per_mass).abs() < 1e-4);
    assert!(
        (grid.total_energy_j() - before - mass * steps as f64 * heat_per_mass).abs()
            < 1e-10 * before.abs()
    );
    let resolved = grid.thermochemical_state([0, 0, 0]).unwrap().unwrap();
    let state = resolved.state();
    assert!((state.density - rho * grid.law().rho_ref_kg_m3).abs() < 1e-12);
    assert!((state.pressure - rho * grid.law().rho_ref_kg_m3 * r * temperature).abs() < 1e-8);
    assert!(state.gamma < reference.state().gamma);
    assert!(state.dynamic_viscosity > reference.state().dynamic_viscosity);
    assert!(state.thermal_conductivity > reference.state().thermal_conductivity);
    assert!(
        (temperature
            - (300.0 + steps as f64 * heat_per_mass / (reference.state().specific_heat_cp - r)))
            .abs()
            > 0.1
    );
}

#[test]
fn hot_source_data_reaches_the_same_native_flow_heat_and_acoustic_state() {
    use fs_material::fluid_transport::FluidTransport;

    let temperature = 2400.0;
    let controls = GasEnergySolve {
        minimum_temperature_k: 2200.0,
        maximum_temperature_k: 2700.0,
        ..controls()
    };
    let gas = ThermochemicalNsfGas::new(
        source::model_with_temperature_limit(-4000.0, 250000.0, 1.0, 3000.0),
        source::point(temperature, 100000.0),
        controls,
    )
    .unwrap();
    let reference = gas.model().evaluate(gas.context()).unwrap();
    let r = reference.state().specific_gas_constant;
    let mut spec = spec(gas);
    spec.temperature_k = temperature;
    spec.temperature_window_k = [
        controls.minimum_temperature_k,
        controls.maximum_temperature_k,
    ];
    let heat_per_mass = 10000.0;
    let source_e = heat_per_mass / (3.0 * r * temperature);
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        spec,
        |_| false,
        |_| (1.0, [0.0; 3], temperature),
        Some(Box::new(move |_, _| (0.0, [0.0; 3], source_e))),
    )
    .unwrap();
    let before = grid.total_energy_j();
    let mass = 32.0 * grid.law().rho_ref_kg_m3 * grid.law().dx_m.powi(3);
    for _ in 0..4 {
        let balance = grid.step().unwrap();
        assert!(balance.residual_rel < 1e-12);
    }
    let resolved = grid.thermochemical_state([0, 0, 0]).unwrap().unwrap();
    let hot = resolved.state();
    hot.validate_coherence().unwrap();
    assert_eq!(
        hot.transport_state().unwrap(),
        resolved.transport_state().unwrap()
    );
    assert!(hot.temperature > temperature);
    assert!(hot.sound_speed > reference.state().sound_speed);
    assert!(hot.gamma < reference.state().gamma);
    assert!(hot.dynamic_viscosity > reference.state().dynamic_viscosity);
    let specific_heat = r
        * (2.0 * (hot.temperature - temperature)
            + 0.0005 * (hot.temperature.powi(2) - temperature.powi(2)));
    assert!((specific_heat - 4.0 * heat_per_mass).abs() < 1e-4);
    assert!((grid.total_energy_j() - before - mass * specific_heat).abs() < 1e-10 * before.abs());
}

#[test]
fn formation_energy_offset_cannot_drive_flow_or_heat_when_density_changes() {
    let mut a =
        NsfGrid3::<D3q27>::new(spec(gas(-4000.0, 250000.0, 1.0)), |_| false, varying).unwrap();
    let mut b =
        NsfGrid3::<D3q27>::new(spec(gas(1_000_000.0, 250000.0, 1.0)), |_| false, varying).unwrap();
    let rho_before = a.cell([0, 0, 0]).0;
    let energy_a = a.total_energy_j();
    let energy_b = b.total_energy_j();
    for _ in 0..5 {
        a.step().unwrap();
        b.step().unwrap();
    }
    assert!((a.cell([0, 0, 0]).0 - rho_before).abs() > 1e-6);
    for x in 0..8 {
        let (ra, ua, ta) = a.cell([x, 0, 0]);
        let (rb, ub, tb) = b.cell([x, 0, 0]);
        assert!((ra - rb).abs() < 1e-9);
        assert!((ta - tb).abs() < 2e-6, "{ta} vs {tb}");
        for axis in 0..3 {
            assert!((ua[axis] - ub[axis]).abs() < 1e-9);
        }
    }
    assert!((a.total_energy_j() - energy_a).abs() < 1e-12 * energy_a.abs());
    assert!((b.total_energy_j() - energy_b).abs() < 1e-12 * energy_b.abs());
    assert!((energy_a - energy_b).abs() < 1e-9 * energy_a.abs());
}

#[test]
fn g3_large_formation_offsets_leave_the_native_thermal_flow_trajectory_unchanged() {
    let mut reference =
        NsfGrid3::<D3q27>::new(spec(gas(0.0, 250000.0, 1.0)), |_| false, varying).unwrap();
    let mut shifted: Vec<_> = [-1.0e12, 1.0e12]
        .into_iter()
        .map(|offset| {
            NsfGrid3::<D3q27>::new(spec(gas(offset, 250000.0, 1.0)), |_| false, varying)
                .unwrap()
        })
        .collect();
    // NASA9 a7 changes the absolute energy convention only. Even the initial
    // sensible energy must agree before density, pressure and heat flux evolve.
    for grid in &shifted {
        assert_eq!(grid.state_hash(), reference.state_hash());
    }
    for _ in 0..5 {
        reference.step().unwrap();
        for grid in &mut shifted {
            grid.step().unwrap();
            assert_eq!(grid.state_hash(), reference.state_hash());
            for x in 0..8 {
                assert_eq!(grid.cell([x, 0, 0]), reference.cell([x, 0, 0]));
            }
        }
    }
}

#[test]
fn g1_small_heat_survives_a_large_material_formation_reference() {
    for heat_per_mass in [-1.0e-3, 1.0e-3] {
        let mut reference_hash = None;
        for offset in [0.0, -1.0e12, 1.0e12] {
            let material = gas(offset, 250000.0, 1.0);
            let r = material.model().evaluate(material.context()).unwrap()
                .state().specific_gas_constant;
            let source_e = heat_per_mass / (3.0 * r * 300.0);
            let mut grid = NsfGrid3::<D3q27>::new_with_sources(
                spec(material), |_| false, |_| (1.0, [0.0; 3], 300.0),
                Some(Box::new(move |_, _| (0.0, [0.0; 3], source_e))),
            ).unwrap();
            for _ in 0..8 { grid.step().unwrap(); }
            let temperature = grid.cell([0, 0, 0]).2;
            // Independent integral of the authored cv/R = 2 + .001 T.
            // Factoring the difference avoids recreating absolute-energy
            // subtraction in the reference calculation itself.
            let heat = r * (temperature - 300.0)
                * (2.0 + 0.0005 * (temperature + 300.0));
            assert!((heat - 8.0 * heat_per_mass).abs() < 1.0e-6, "{offset}: {heat}");
            let expected = reference_hash.get_or_insert(grid.state_hash());
            assert_eq!(grid.state_hash(), *expected);
        }
    }
}

#[test]
fn conductivity_data_changes_native_energy_exchange_and_worker_count_does_not() {
    let initial = |x: [usize; 3]| {
        let t = 300.0 + 2.0 * (2.0 * core::f64::consts::PI * x[0] as f64 / 8.0).cos();
        (300.0 / t, [0.0; 3], t) // uniform initial physical pressure
    };
    let mut one =
        NsfGrid3::<D3q27>::new(spec(gas(-4000.0, 250000.0, 1.0)), |_| false, initial).unwrap();
    let mut parallel_spec = spec(gas(-4000.0, 250000.0, 1.0));
    parallel_spec.workers = 3;
    let mut parallel = NsfGrid3::<D3q27>::new(parallel_spec, |_| false, initial).unwrap();
    let mut conducting =
        NsfGrid3::<D3q27>::new(spec(gas(-4000.0, 250000.0, 4.0)), |_| false, initial).unwrap();
    for _ in 0..5 {
        one.step().unwrap();
        parallel.step().unwrap();
        conducting.step().unwrap();
    }
    assert_eq!(one.state_hash(), parallel.state_hash());
    assert!((one.cell([0, 0, 0]).2 - conducting.cell([0, 0, 0]).2).abs() > 1e-4);
    assert!(conducting.cell([0, 0, 0]).2 < one.cell([0, 0, 0]).2);
}

#[test]
fn actual_local_pressure_support_refuses_initial_and_future_states_atomically() {
    assert!(matches!(
        NsfGrid3::<D3q27>::new(
            spec(gas(-4000.0, 110000.0, 1.0)),
            |_| false,
            |_| (1.2, [0.0; 3], 300.0)
        ),
        Err(NsfError::Gas(_))
    ));
    let source = gas(-4000.0, 110000.0, 1.0);
    let reference = source.model().evaluate(source.context()).unwrap();
    let r = reference.state().specific_gas_constant;
    let source_e = 100000.0 / (3.0 * r * 300.0);
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        spec(source),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(move |_, _| (0.0, [0.0; 3], source_e))),
    )
    .unwrap();
    let before = (
        grid.state_hash(),
        grid.total_energy_j(),
        grid.cell([0, 0, 0]),
        grid.cumulative_ledger(),
    );
    assert!(matches!(grid.step(), Err(NsfError::Gas(_))));
    assert_eq!(grid.steps(), 0);
    assert_eq!(
        before,
        (
            grid.state_hash(),
            grid.total_energy_j(),
            grid.cell([0, 0, 0]),
            grid.cumulative_ledger()
        )
    );
    assert!(grid.step().is_err());
    assert_eq!(before.0, grid.state_hash());
}

#[test]
fn local_temperature_dependent_gamma_controls_mach_admission() {
    let source = gas(-4000.0, 250000.0, 1.0);
    let gamma_reference = source
        .model()
        .evaluate(source.context())
        .unwrap()
        .state()
        .gamma;
    let speed = 0.1999 * (gamma_reference * 550.0 / 900.0).sqrt();
    assert!(matches!(
        NsfGrid3::<D3q27>::new(spec(source), |_| false, |_| (1.0, [speed, 0.0, 0.0], 550.0)),
        Err(NsfError::Mach { .. })
    ));
}

#[test]
fn perfect_gas_ingress_retains_legacy_arithmetic_and_observations() {
    let gas = GasSpec::dry_air_ussa1976();
    let mut direct = NsfGrid3::<D3q27>::new(spec(gas.clone()), |_| false, varying).unwrap();
    let mut explicit = NsfGrid3::<D3q27>::new(
        spec(NsfGasMedium::CaloricallyPerfect(gas)),
        |_| false,
        varying,
    )
    .unwrap();
    assert!(direct.thermochemical_state([0, 0, 0]).unwrap().is_none());
    for _ in 0..5 {
        assert_eq!(direct.step().unwrap(), explicit.step().unwrap());
        assert_eq!(direct.state_hash(), explicit.state_hash());
    }
}

#[test]
fn signed_heat_sources_close_native_energy_without_becoming_viscous_heat() {
    for thermochemical in [false, true] {
        let mut expected_replay = None;
        for workers in [1, 3] {
            let medium = if thermochemical {
                NsfGasMedium::Thermochemical(gas(-4000.0, 250000.0, 1.0))
            } else {
                NsfGasMedium::CaloricallyPerfect(GasSpec::dry_air_ussa1976())
            };
            let mut input = spec(medium);
            input.workers = workers;
            let mut grid = NsfGrid3::<D3q27>::new_with_sources(
                input,
                |x| x[2] == 1,
                |_| (1.0, [0.0; 3], 300.0),
                Some(Box::new(|x, time| {
                    assert_eq!(x[2], 0, "solid cells must not request a source");
                    (0.0, [0.0; 3], 0.002 * (2.0 - time))
                })),
            )
            .unwrap();
            let initial = grid.total_energy_j();
            let unit = grid.law().energy_unit_j();
            let mut supplied = 0.0;
            for step in 0..4 {
                let before = grid.total_energy_j();
                let ledger = grid.step().unwrap();
                // Sixteen fluid cells; trapezoidal integration is exact for
                // this independently prescribed affine power history.
                let expected = 16.0 * unit * 0.002 * (1.5 - step as f64);
                supplied += expected;
                assert!((ledger.source_energy_j - expected).abs() < 1e-13 * initial);
                assert!((grid.total_energy_j() - before - expected).abs() < 1e-12 * initial);
                assert!(ledger.residual_rel < 1e-12);
                assert!(ledger.viscous_heating_j.abs() < 1e-12 * initial);
                assert_eq!(ledger.unresolved_internal_exchange_j, 0.0);
                assert!(
                    (grid.cumulative_ledger().source_energy_j - supplied).abs() < 1e-12 * initial
                );
            }
            assert!((grid.total_energy_j() - initial).abs() < 1e-12 * initial);
            let replay = (grid.state_hash(), grid.cumulative_ledger());
            if let Some(expected) = expected_replay {
                assert_eq!(replay, expected);
            } else {
                expected_replay = Some(replay);
            }
        }
    }
}

#[test]
fn mechanical_source_work_is_never_relabelled_as_viscous_heat() {
    let force: f64 = 0.002;
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        spec(GasSpec::dry_air_ussa1976()),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(move |_, time| {
            // Uniform acceleration: u = force*time, hence supplied power
            // is force*u and all supplied energy becomes kinetic energy.
            (0.0, [force, 0.0, 0.0], force * force * time)
        })),
    )
    .unwrap();
    let initial = grid.total_energy_j();
    for step in 0..4 {
        let before_kinetic = grid.kinetic_energy_j();
        let ledger = grid.step().unwrap();
        let expected = 32.0 * grid.law().energy_unit_j() * force * force * (step as f64 + 0.5);
        assert!((ledger.source_energy_j - expected).abs() < 1e-13 * initial);
        assert!((grid.kinetic_energy_j() - before_kinetic - expected).abs() < 1e-12 * initial);
        assert_eq!(ledger.viscous_heating_j, 0.0);
        assert!((ledger.unresolved_internal_exchange_j + expected).abs() < 1e-12 * initial);
        assert!(ledger.residual_rel < 1e-12);
        assert!((grid.cell([0, 0, 0]).2 - 300.0).abs() < 1e-8);
    }
}

#[test]
fn mass_source_energy_is_a_supply_not_dissipation() {
    let gas = GasSpec::dry_air_ussa1976();
    let velocity: f64 = 0.02;
    let mass: f64 = 0.001;
    let energy_per_mass = (1.0 / 3.0) / (gas.gamma - 1.0) + 0.5 * velocity * velocity;
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        spec(gas),
        |_| false,
        |_| (1.0, [velocity, 0.0, 0.0], 300.0),
        Some(Box::new(move |_, _| {
            (mass, [mass * velocity, 0.0, 0.0], mass * energy_per_mass)
        })),
    )
    .unwrap();
    let initial = grid.total_energy_j();
    let ledger = grid.step().unwrap();
    let unit = 32.0 * grid.law().energy_unit_j();
    assert!((ledger.source_energy_j - unit * mass * energy_per_mass).abs() < 1e-13 * initial);
    assert!(ledger.residual_rel < 1e-12);
    assert_eq!(ledger.viscous_heating_j, 0.0);
    assert!(
        (ledger.unresolved_internal_exchange_j + unit * mass * 0.5 * velocity * velocity).abs()
            < 1e-12 * initial
    );
    let (density, actual_velocity, temperature) = grid.cell([0, 0, 0]);
    assert!((density - 1.0 - mass).abs() < 1e-12);
    assert!((actual_velocity[0] - velocity).abs() < 1e-12);
    assert!((temperature - 300.0).abs() < 1e-8);
}

fn named_gas() -> NamedMixtureNsfGas {
    NamedMixtureNsfGas::new(
        vec![GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::CarbonDioxide,
            GasSpecies::WaterVapor],
        Composition::new(CompositionBasis::MoleFraction, vec![0.755, 0.21, 0.025, 0.01]).unwrap(),
        300.0,
        100000.0,
    ).unwrap()
}

fn named_spec(gas: NamedMixtureNsfGas) -> NsfSpec<NamedMixtureNsfGas> {
    let mut input = spec(gas);
    input.temperature_window_k = [285.0, 320.0];
    input
}

#[test]
fn g0_named_mixture_keeps_species_mass_basis_and_original_local_material() {
    let species = vec![GasSpecies::Helium, GasSpecies::Oxygen, GasSpecies::CarbonDioxide];
    let composition = Composition::new(CompositionBasis::MassFraction, vec![0.4, 0.5, 0.1]).unwrap();
    let source = NamedMixtureNsfGas::new(species.clone(), composition.clone(), 300.0, 100000.0).unwrap();
    let reference = resolve_gas_mixture(&species, &composition, 300.0, 100000.0).unwrap();
    assert_eq!(source.reference_state(), &reference);
    let grid = NsfGrid3::<D3q27>::new(named_spec(source), |x| x[2] == 1, varying).unwrap();
    let NsfGasMedium::NamedMixture(retained) = grid.medium() else {
        panic!("the grid must retain the original named-mixture owner");
    };
    assert_eq!(retained.species(), species.as_slice());
    assert_eq!(retained.composition(), &composition);
    assert_eq!(grid.law().gas_constant, reference.specific_gas_constant);
    assert_eq!(grid.law().rho_ref_kg_m3, reference.density);
    for x in 0..8 {
        let state = grid.gas_state([x, 0, 0]).unwrap();
        let expected = resolve_gas_mixture(&species, &composition, state.temperature, state.pressure).unwrap();
        assert_eq!(state, expected);
        let density = grid.cell([x, 0, 0]).0 * grid.law().rho_ref_kg_m3;
        assert!((state.density / density - 1.0).abs() < 1e-14);
        assert!(grid.thermochemical_state([x, 0, 0]).unwrap().is_none());
    }
    assert!(grid.gas_state([0, 0, 1]).is_err());
    assert!(grid.gas_state([8, 0, 0]).is_err());
}

#[test]
fn g1_g5_named_mixture_heating_closes_source_cv_energy_and_replays_across_workers() {
    let reference = *named_gas().reference_state();
    let cv = reference.specific_heat_cp - reference.specific_gas_constant;
    let heat_per_mass = 2.0 * cv;
    let source_energy = heat_per_mass / (3.0 * reference.specific_gas_constant * reference.temperature);
    let mut replay = None;
    for (workers, transport_points) in [(1, 2), (3, 17)] {
        let mut input = named_spec(named_gas());
        input.workers = workers;
        input.transport_points = transport_points;
        let mut grid = NsfGrid3::<D3q27>::new_with_sources(
            input,
            |_| false,
            |_| (1.0, [0.0; 3], 300.0),
            Some(Box::new(move |_, _| (0.0, [0.0; 3], source_energy))),
        ).unwrap();
        let initial = grid.total_energy_j();
        let mass = 32.0 * reference.density * grid.law().dx_m.powi(3);
        for _ in 0..4 {
            assert!(grid.step().unwrap().residual_rel < 1e-12);
        }
        let state = grid.gas_state([0, 0, 0]).unwrap();
        assert!((state.temperature - 308.0).abs() < 1e-9);
        assert!((grid.total_energy_j() - initial - mass * 4.0 * heat_per_mass).abs() < 1e-12 * initial);
        assert!((state.pressure / reference.pressure - 308.0 / 300.0).abs() < 1e-12);
        assert_eq!(state.specific_heat_cp, reference.specific_heat_cp);
        assert_eq!(state.specific_gas_constant, reference.specific_gas_constant);
        assert!(state.dynamic_viscosity > reference.dynamic_viscosity);
        assert!(state.thermal_conductivity > reference.thermal_conductivity);
        assert!(state.sound_speed > reference.sound_speed);
        let result = (grid.state_hash(), grid.cumulative_ledger());
        if let Some(expected) = replay { assert_eq!(result, expected); }
        replay = Some(result);
    }
}

#[test]
fn g1_g5_named_mixture_local_transport_and_sgs_heat_ignore_diagnostic_table_spacing() {
    let mut grids = Vec::new();
    for (workers, transport_points) in [(1, 2), (3, 17)] {
        let mut input = named_spec(named_gas());
        input.workers = workers;
        input.transport_points = transport_points;
        grids.push(NsfGrid3::<D3q27>::new(input, |_| false, |x| {
            let phase = core::f64::consts::TAU * (x[0] as f64 + 0.5) / 8.0;
            let temperature = 300.0 + 3.0 * phase.cos();
            (300.0 / temperature, [0.0, 0.02 * phase.sin(), 0.0], temperature)
        }).unwrap().with_sgs(NsfSgsConfig {
            closure: SgsConfig {
                model: SgsModel::Smagorinsky { cs: 0.17 },
                filter_width_lu: 1.0,
                maximum_eddy_to_molecular: 10.0,
                usage: SgsUse::Exploration,
            },
            turbulent_prandtl: 0.85,
        }).unwrap());
    }
    let initial = grids[0].total_energy_j();
    for _ in 0..4 {
        assert_eq!(grids[0].state_hash(), grids[1].state_hash());
        for x in 0..8 {
            let location = [x, 0, 0];
            let state = grids[0].gas_state(location).unwrap();
            let density = grids[0].cell(location).0;
            let transport = grids[0].sgs_transport(location).unwrap().unwrap();
            assert!(transport.eddy_viscosity_lu > 0.0);
            let expected = density * state.specific_heat_cp / state.specific_gas_constant
                * transport.eddy_viscosity_lu / 0.85;
            assert!((transport.turbulent_conductivity_lu / expected - 1.0).abs() < 1e-12);
        }
        let a = grids[0].step().unwrap();
        let b = grids[1].step().unwrap();
        assert_eq!(a, b);
        assert!(a.residual_rel < 1e-12);
        assert!((grids[0].total_energy_j() / initial - 1.0).abs() < 2e-12);
    }
    assert_eq!(grids[0].state_hash(), grids[1].state_hash());
}

#[test]
fn g1_named_mixture_warm_open_flow_uses_the_same_sensible_energy() {
    let mut input = named_spec(named_gas());
    input.periodic[0] = false;
    let mut grid = NsfGrid3::<D3q27>::new(
        input, |_| false, |_| (300.0 / 305.0, [0.02, 0.0, 0.0], 305.0),
    ).unwrap();
    let exterior = NsfOpenState {
        density_kg_m3: grid.law().rho_ref_kg_m3 * 300.0 / 305.0,
        velocity_m_s: [0.02 * grid.law().dx_m / grid.law().dt_s, 0.0, 0.0],
        temperature_k: 305.0,
    };
    for high in [false, true] { grid.set_open_face(0, high, Some(exterior)).unwrap(); }
    for _ in 0..8 { assert!(grid.step().unwrap().residual_rel < 1e-12); }
    for x in 0..8 {
        let (density, velocity, temperature) = grid.cell([x, 0, 0]);
        assert!((density - 300.0 / 305.0).abs() < 1e-10);
        assert!((velocity[0] - 0.02).abs() < 1e-10);
        assert!((temperature - 305.0).abs() < 1e-8);
    }
}

#[test]
fn g4_named_mixture_domain_and_condensation_refuse_without_publishing() {
    let gas = named_gas();
    for (temperature, pressure) in [(270.0, 100000.0), (330.0, 100000.0),
        (300.0, 9000.0), (300.0, 210000.0)] {
        assert!(matches!(gas.at(temperature, pressure), Err(NsfError::Gas(_))));
    }
    let mut input = named_spec(gas.clone());
    input.temperature_k = 301.0;
    assert!(matches!(NsfGrid3::<D3q27>::new(input, |_| false, |_| (1.0, [0.0; 3], 301.0)),
        Err(NsfError::Input(_))));
    let mut input = named_spec(gas.clone());
    input.temperature_window_k[0] = 270.0;
    assert!(matches!(NsfGrid3::<D3q27>::new(input, |_| false, |_| (1.0, [0.0; 3], 300.0)),
        Err(NsfError::Gas(_))));
    assert!(matches!(NsfGrid3::<D3q27>::new(named_spec(gas.clone()), |_| false,
        |_| (2.1, [0.0; 3], 300.0)), Err(NsfError::Gas(_))));

    let reference = gas.reference_state();
    let excessive_heat = 30.0 * (reference.specific_heat_cp - reference.specific_gas_constant)
        / (3.0 * reference.specific_gas_constant * reference.temperature);
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(named_spec(gas), |_| false,
        |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(move |_, _| (0.0, [0.0; 3], excessive_heat))),
    ).unwrap();
    let before = (grid.state_hash(), grid.cumulative_ledger(), grid.steps());
    assert!(matches!(grid.step(), Err(NsfError::Gas(_))));
    assert_eq!((grid.state_hash(), grid.cumulative_ledger(), grid.steps()), before);

    let wet = NamedMixtureNsfGas::new(
        vec![GasSpecies::Nitrogen, GasSpecies::WaterVapor],
        Composition::new(CompositionBasis::MoleFraction, vec![0.97, 0.03]).unwrap(),
        300.0, 100000.0,
    ).unwrap();
    let mut input = named_spec(wet);
    input.temperature_window_k = [299.0, 310.0];
    input.periodic[0] = false;
    let mut grid = NsfGrid3::<D3q27>::new(input, |_| false, |_| (1.0, [0.0; 3], 300.0)).unwrap();
    let before = grid.state_hash();
    let compressed = NsfOpenState {
        density_kg_m3: 1.3 * grid.law().rho_ref_kg_m3,
        velocity_m_s: [0.0; 3], temperature_k: 300.0,
    };
    assert!(matches!(grid.set_open_face(0, false, Some(compressed)), Err(NsfError::Gas(_))));
    assert_eq!(grid.state_hash(), before);
}
