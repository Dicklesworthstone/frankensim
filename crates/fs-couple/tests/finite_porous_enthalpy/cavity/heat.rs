//! Physical heat exchange on retained finite compression histories.
#[path = "heat/real_gas.rs"]
mod real_gas;

#[path = "heat/source.rs"]
mod source;

use super::*;
use fs_couple::solid_enthalpy::finite_plastic::porous::cavity::heat::{
    CavityHeatError, CavityHeatObserverError, CavityHeatSolve, CavityWallHeatLink,
    ThermalPorousGasCavity, ThermalPorousGasSystem,
};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_solver::NewtonKrylovConfig;

/// Original source fixture with explicitly zero thermoelastic/softening terms.
/// This isolates an independent two-capacity heat oracle, without replacing
/// either the native gas calorics or the porous solid/thermal owners.
fn thermally_constant_fixture() -> Fixture {
    let mut constitutive = model();
    for key in [
        "shear_modulus_slope_pa_k", "bulk_modulus_slope_pa_k",
        "yield_stress_slope_pa_k", "hardening_modulus_slope_pa_k",
        "voce_saturation_slope_pa_k", "logarithmic_expansion_per_k",
        "logarithmic_expansion_slope_per_k2",
    ] {
        constitutive.parameters.get_mut(key).unwrap().value = 0.0;
    }
    let mut claims = ClaimSet::new();
    claims.insert_claim(scalar(CP_NAME, CP_DIMS, CP)).unwrap();
    claims.insert_claim(scalar(DENSITY_PROPERTY, Density::DIMS, RHO)).unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored constant thermal porous heat exchange".into(),
            phase: "solid".into(), process: "numerical regression".into(), revision: 0,
        },
        claims, vec![constitutive],
    ).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "temperature",
            lower: &source_context(280.0), upper: &source_context(420.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid, -20.0 * CP, 1e-6,
    ).unwrap();
    Fixture { card, caloric }
}

fn solve_heat() -> CavityHeatSolve {
    CavityHeatSolve {
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0, maximum_temperature_k: 600.0,
            energy_tolerance_j_kg: 1e-7, maximum_iterations: 128,
        },
        maximum_trials: 128, maximum_newton_iterations: 20, difference_step_j: 1e-5,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9, relative_tolerance: 0.0, ..Default::default()
        },
        energy_tolerance_j: 1e-7,
    }
}

fn cavity_entry(cavity: Arc<ClosedTetCavity>, gas: ClosedGasParcel) -> ThermalPorousGasCavity {
    ThermalPorousGasCavity {
        cavity, gas, exterior_pressure_pa: 2e5,
        entropy: GasEntropySolve {
            minimum_temperature_k: 200.0, maximum_temperature_k: 600.0,
            entropy_tolerance_j_kg_k: 1e-10, maximum_iterations: 128,
        },
    }
}

/// Two distinct declared patches per cell, with total G/C = rate.
fn links(thermal: &DepositionSimulation<'_>, rate: f64) -> Vec<CavityWallHeatLink> {
    thermal.cells().values().enumerate().flat_map(|(i, cell)| {
        [0.4, 0.6].into_iter().enumerate().map(move |(patch, fraction)| {
            CavityWallHeatLink {
                id: (2 * i + patch) as u64, cavity: 0, cell: cell.id,
                conductance_w_k: fraction * rate * cell.mass_kg * CP,
            }
        })
    }).collect()
}

#[derive(Clone, Default)]
struct Times(Vec<f64>);
impl DepositionObserver for Times {
    type Error = DepositionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.0.clear();
        Ok(())
    }
    fn step(
        &mut self, _: f64, dt: f64, _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        self.0.push(dt);
        Ok(())
    }
}

#[test]
fn rigid_finite_gas_wall_heat_matches_independent_capacity_balance() {
    context(|cx| {
        let fixture = thermally_constant_fixture();
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let cold_gas = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.0, 1e24);
        let (gas, _) = cold_gas.at_temperature(360.0).unwrap();
        let original_gas = gas.clone();
        let rate = 0.4;
        let system = ThermalPorousGasSystem::try_new(
            solid, &thermal, vec![cavity_entry(Arc::new(cavity), gas)],
            links(&thermal, rate), cx,
        ).unwrap();
        let observer = Times::default();
        let external = 0.02;
        let next = system.advance(
            &thermal, &observer, 1.0, &affine(&nodes, 1.0), &[], [0.0; 3], None,
            cavity_settings(), controls(), thermal_config(), config(),
            CouplingMethod::RelaxedPicard, &[external], &[0.0], solve_heat(), cx,
        ).unwrap();
        assert!(next.observer.0.len() > 4);
        close(next.observer.0.iter().sum(), 1.0, 1e-15);
        // Independent closed-form explicit wall response on the actual native
        // time grid, coupled to one implicit finite-gas energy equation.
        let decay: f64 = next.observer.0.iter().map(|dt| 1.0 - rate * dt).product();
        let wall_capacity: f64 = thermal.cells().values().map(|cell| cell.mass_kg * CP).sum();
        let gas_capacity = original_gas.mass_kg()
            * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028 * 2.5;
        let exchange_capacity = wall_capacity * (1.0 - decay);
        let expected_gas = (gas_capacity * 360.0 + exchange_capacity * 300.0 + external)
            / (gas_capacity + exchange_capacity);
        let expected_wall = 300.0 + (expected_gas - 300.0) * (1.0 - decay);
        let gas = &next.system.cavities()[0].gas;
        close(gas.state().state().temperature, expected_gas, 3e-7);
        close(next.balance.gases[0].heat_stage_temperature_k, expected_gas, 3e-7);
        close(next.balance.gases[0].adiabatic_energy_change_j, 0.0, 0.0);
        close(next.balance.total_gas_energy_change_j,
            gas_capacity * (expected_gas - 360.0), 3e-9);
        close(next.balance.gas_to_solid_heat_j,
            wall_capacity * (expected_wall - 300.0), 3e-9);
        close(next.balance.external_solid_heat_j, 0.0, 1e-9);
        close(next.balance.native.applied_heat_j, 0.0, 1e-9);
        close(next.balance.total_gas_energy_change_j
            + next.balance.native.thermal_audit.stored_energy_change_j, external, 3e-9);
        close(next.balance.combined_residual_j, 0.0, solve_heat().energy_tolerance_j);
        assert!(next.balance.interface_heat_defect_j < 1e-12);
        assert_eq!(next.balance.links.len(), 2 * tets.len());
        for (i, cell) in thermal.cells().values().enumerate() {
            close(next.thermal.temperature_k(cell.id).unwrap(), expected_wall, 3e-7);
            close(next.balance.links[2 * i].heat_into_solid_j,
                0.4 * cell.mass_kg * CP * (expected_wall - 300.0), 1e-9);
            close(next.balance.links[2 * i + 1].heat_into_solid_j,
                0.6 * cell.mass_kg * CP * (expected_wall - 300.0), 1e-9);
        }
        assert_eq!(gas.volume_m3().to_bits(), original_gas.volume_m3().to_bits());
        assert_eq!(gas.mass_kg().to_bits(), original_gas.mass_kg().to_bits());
        assert_eq!(gas.model(), original_gas.model());
        assert_eq!(system.cavities()[0].gas, original_gas);
        assert!(observer.0.is_empty());
        assert_eq!(thermal.time_s(), 0.0);
        assert!(thermal.cells().values().all(|cell| cell.specific_h_j_kg == 0.0));
    });
}

/// Independent variable-Cv isentrope from THIS interval's heated actual state.
fn epoch_temperature(initial: f64, volume_ratio: f64, slope: f64) -> f64 {
    let (mut lo, mut hi) = (200.0, 600.0);
    for _ in 0..100 {
        let t = 0.5 * (lo + hi);
        if 2.5 * (t / initial).ln() + slope * (t - initial) + volume_ratio.ln() < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn heated_compression_retains_gas_compaction_inertia_and_atomic_failures() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx)
            .with_inertia(&vec![[0.001, -0.0005, 0.0003]; nodes.len()], cx).unwrap();
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let cold_gas = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.001, 1e24);
        let (gas, _) = cold_gas.at_temperature(360.0).unwrap();
        let initial_gas = gas.clone();
        let system = ThermalPorousGasSystem::try_new(
            solid, &thermal, vec![cavity_entry(Arc::new(cavity), gas)],
            links(&thermal, 0.005), cx,
        ).unwrap();
        let first = system.advance(
            &thermal, &(), 1.0, &affine(&nodes, 0.94), &[], [0.0; 3], None,
            cavity_settings(), controls(), thermal_config(), config(),
            CouplingMethod::IqnIls(IqnIlsConfig::default()),
            &[0.0], &[0.0], solve_heat(), cx,
        ).unwrap();
        assert!(first.system.solid().mechanics().element_states().iter().all(|state| state.alpha > 0.0));
        assert!(first.balance.native.inertia.is_some());
        assert!(first.balance.gas_to_solid_heat_j > 0.0);
        assert!(first.balance.native.applied_heat_j > 0.0);
        let before_volume = first.system.cavities()[0].gas.volume_m3();
        close(first.system.cavities()[0].gas.state().state().temperature,
            epoch_temperature(first.balance.gases[0].heat_stage_temperature_k,
                before_volume / initial_gas.volume_m3(), 0.001), 3e-7);
        let old_cells = first.thermal.cells().clone();
        let old_gas = first.system.cavities()[0].gas.clone();
        let old_states = first.system.solid().mechanics().element_states().to_vec();
        let old_displacement = first.system.solid().mechanics().displacement_m().to_vec();
        let old_velocity = first.system.solid().mechanics().velocity_mps().unwrap().to_vec();
        let supports = affine(&nodes, 0.90);
        let refusal = LateRefusal { completed_substeps: 0 };
        let refused = first.system.advance(
            &first.thermal, &refusal, 2.0, &supports, &[], [0.0; 3], None,
            cavity_settings(), controls(), thermal_config(), config(),
            CouplingMethod::RelaxedPicard, &[0.0], &[0.0], solve_heat(), cx,
        );
        assert!(matches!(refused, Err(CavityHeatError::Feedback(
            FiniteThermoplasticError::Observer(CavityHeatObserverError::Observer(
                DepositionError::Cancelled
            ))
        ))));
        let limited = first.system.advance(
            &first.thermal, &(), 2.0, &supports, &[], [0.0; 3], None,
            cavity_settings(), controls(), thermal_config(), config(),
            CouplingMethod::RelaxedPicard, &[0.0], &[0.0],
            CavityHeatSolve { maximum_trials: 1, ..solve_heat() }, cx,
        );
        assert!(matches!(limited, Err(CavityHeatError::TrialBudget)));
        let unsupported = first.system.advance(
            &first.thermal, &(), 2.0, &supports, &[], [0.0; 3], None,
            cavity_settings(), controls(), thermal_config(), config(),
            CouplingMethod::RelaxedPicard, &[0.0], &[1e12], solve_heat(), cx,
        );
        assert!(matches!(unsupported, Err(CavityHeatError::Gas(_))));
        let second = first.system.advance(
            &first.thermal, &(), 2.0, &supports, &[], [0.0; 3], None,
            cavity_settings(), controls(), thermal_config(), config(),
            CouplingMethod::RelaxedPicard, &[0.0], &[0.0], solve_heat(), cx,
        ).unwrap();
        let final_gas = &second.system.cavities()[0].gas;
        let final_volume = second.system.cavities()[0].cavity.evaluate(
            second.system.solid().mechanics().displacement_m().as_flattened(), cx,
        ).unwrap().volume_m3();
        assert_eq!(final_gas.volume_m3().to_bits(), final_volume.to_bits());
        close(final_gas.state().state().temperature,
            epoch_temperature(second.balance.gases[0].heat_stage_temperature_k,
                final_volume / before_volume, 0.001), 3e-7);
        close(second.balance.total_gas_energy_change_j,
            gas_energy(&old_gas, 0.001, old_gas.state().state().temperature,
                final_gas.state().state().temperature), 3e-8);
        close(second.balance.external_solid_heat_j, 0.0, 2e-8);
        close(second.balance.combined_residual_j, 0.0, solve_heat().energy_tolerance_j);
        assert_eq!(final_gas.mass_kg().to_bits(), initial_gas.mass_kg().to_bits());
        assert_eq!(final_gas.model(), initial_gas.model());
        for (i, state) in second.system.solid().mechanics().element_states().iter().enumerate() {
            assert!(state.alpha >= old_states[i].alpha);
            close(second.thermal.cells()[&CellId(100 + i as u64)].mass_kg,
                RHO * volume(&nodes, tets[i]), 0.0);
            assert!(std::ptr::eq(
                second.system.solid().mechanics().families()[i].caloric_source(), &fixture.caloric,
            ));
        }
        assert_eq!(refusal.completed_substeps, 0);
        assert_eq!(first.thermal.cells(), &old_cells);
        assert_eq!(first.system.cavities()[0].gas, old_gas);
        assert_eq!(first.system.solid().mechanics().element_states(), old_states);
        assert_eq!(first.system.solid().mechanics().displacement_m(), old_displacement);
        assert_eq!(first.system.solid().mechanics().velocity_mps().unwrap(), old_velocity);
        assert_eq!(first.system.time_s(), 1.0);
        assert_eq!(second.system.time_s(), 2.0);
    });
}

#[test]
fn supplied_gas_heat_changes_an_actually_free_porous_wall_through_native_pressure() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
        let (surface, cavity) = boundaries(&nodes, &tets, cx);
        let stretch: f64 = 0.94;
        let temperature = adiabatic_temperature(stretch, CP);
        let pressure = -hydrostatic(stretch, temperature).first_piola_pa / stretch.powi(2);
        let exterior = 2e5;
        let gas = gas_fixture::gas(
            cavity.reference_volume_m3(), (pressure + exterior) * stretch.powf(4.2),
            0.0, 1e24,
        );
        let fields: Vec<_> = (0..surface.faces().len()).map(|face| {
            if cavity.solid_face_indices().contains(&face) { 0.0 } else { pressure }
        }).collect();
        let supports: Vec<_> = affine(&nodes, stretch).into_iter()
            .filter(|(dof, _)| *dof != 15).collect();
        let system = ThermalPorousGasSystem::try_new(
            solid, &thermal, vec![cavity_entry(Arc::new(cavity), gas)],
            links(&thermal, 0.005), cx,
        ).unwrap();
        let run = |external: f64| system.advance(
            &thermal, &(), 1.0, &supports, &[], [0.0; 3],
            Some((&surface, TetPressure::PerFace(&fields))),
            cavity_settings(), controls(), thermal_config(), config(),
            CouplingMethod::RelaxedPicard, &[external], &[0.0], solve_heat(), cx,
        ).unwrap();
        let unheated = run(0.0);
        let heated = run(0.005);
        assert!(heated.balance.native.equilibrium.solid.steps[0].linear_iterations > 0);
        assert!(heated.system.solid().mechanics().displacement_m()[5][0]
            > unheated.system.solid().mechanics().displacement_m()[5][0] + 1e-10);
        assert!(heated.system.cavities()[0].gas.volume_m3()
            > unheated.system.cavities()[0].gas.volume_m3());
        for next in [&unheated, &heated] {
            close(next.balance.gases[0].heat_residual_j, 0.0, solve_heat().energy_tolerance_j);
            close(next.balance.combined_residual_j, 0.0, solve_heat().energy_tolerance_j);
            assert!(next.balance.native.material_energy.pressure.is_some());
        }
        assert_eq!(system.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
