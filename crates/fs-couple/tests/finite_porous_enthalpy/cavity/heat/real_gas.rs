//! Dense-gas heat/storage/acoustics on the same retained porous wall owner.
//! Synthetic PR critical data and calorics are numerical regression sources,
//! not a measured methane material card.
use super::*;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    energy::parcel::ClosedRealGasParcel,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

const MOLAR_MASS: f64 = 0.028;
const TC: f64 = 190.564;
const PC: f64 = 4.5992e6;
const OMEGA: f64 = 0.01142;
const SLOPE: f64 = 0.001;
const P0: f64 = 5e6;

fn real_parcel(volume: f64, datum: f64) -> ClosedRealGasParcel {
    let original = super::super::gas_fixture::gas(volume, P0, SLOPE, datum);
    let model = ThermochemicalPengRobinsonModel::new(
        original.model().clone(),
        vec![PengRobinsonSpeciesBinding {
            species: SpeciesId::new("synthetic-cavity-gas").unwrap(),
            critical_temperature: Temperature::new(TC),
            critical_pressure: Pressure::new(PC),
            acentric_factor: OMEGA,
        }],
        vec![],
    ).unwrap();
    ClosedRealGasParcel::new(
        model, original.point().clone(), volume, PengRobinsonPhaseV1::Vapor,
    ).unwrap()
}

#[derive(Clone, Copy)]
struct Pr {
    pressure: f64,
    residual_u: f64,
    attractive_s: f64,
    bulk: f64,
}

fn b() -> f64 { 0.07780 * R * TC / PC }

/// Independent pure-fluid Helmholtz derivatives at original molar volume.
/// No native EOS, heat inversion, entropy solve or acoustic owner is called.
fn pr(t: f64, v: f64) -> Pr {
    let k = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA * OMEGA;
    let a0 = 0.45724 * R * R * TC * TC / PC;
    let root = 1.0 + k * (1.0 - (t / TC).sqrt());
    let a = a0 * root * root;
    let da = -a0 * k * root / (TC * t).sqrt();
    let dda = a0 * k * (1.0 + k) / (2.0 * TC.sqrt() * t.powf(1.5));
    let b = b();
    let root_two = 2.0_f64.sqrt();
    let integral = ((v + (1.0 + root_two) * b) / (v + (1.0 - root_two) * b)).ln()
        / (2.0 * root_two * b);
    let denominator = v * v + 2.0 * b * v - b * b;
    let pressure = R * t / (v - b) - a / denominator;
    let p_t = R / (v - b) - da / denominator;
    let p_v = -R * t / (v - b).powi(2) + 2.0 * a * (v + b) / denominator.powi(2);
    let cv = R * (2.5 + SLOPE * t) + t * dda * integral;
    Pr {
        pressure, residual_u: (t * da - a) * integral, attractive_s: da * integral,
        bulk: -v * (p_v - t * p_t * p_t / cv),
    }
}

fn initial_molar_volume() -> f64 {
    let (mut lo, mut hi) = (1.001 * b(), 20.0 * R * 300.0 / P0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if pr(300.0, mid).pressure > P0 { lo = mid; } else { hi = mid; }
    }
    0.5 * (lo + hi)
}

fn energy(t0: f64, v0: f64, t1: f64, v1: f64) -> f64 {
    R * (t1 - t0) * (2.5 + 0.5 * SLOPE * (t1 + t0))
        + pr(t1, v1).residual_u - pr(t0, v0).residual_u
}

fn entropy(t0: f64, v0: f64, t1: f64, v1: f64) -> f64 {
    R * (2.5 * (t1 / t0).ln() + SLOPE * (t1 - t0) + ((v1 - b()) / (v0 - b())).ln())
        + pr(t1, v1).attractive_s - pr(t0, v0).attractive_s
}

fn entry(cavity: Arc<ClosedTetCavity>, gas: ClosedRealGasParcel)
    -> ThermalPorousGasCavity<ClosedRealGasParcel>
{
    ThermalPorousGasCavity {
        cavity, gas, exterior_pressure_pa: 2e5,
        entropy: GasEntropySolve {
            minimum_temperature_k: 200.0, maximum_temperature_k: 600.0,
            entropy_tolerance_j_kg_k: 1e-10, maximum_iterations: 128,
        },
    }
}

#[test]
fn dense_gas_wall_heat_uses_original_residual_calorics_and_true_acoustic_bulk() {
    context(|cx| {
        let fixture = thermally_constant_fixture();
        let (nodes, tets) = shell();
        let mut without_datum: Option<[f64; 4]> = None;
        for datum in [0.0, 1e24] {
            let thermal = original_thermal(&fixture, &nodes, &tets);
            let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
            let (_, cavity) = boundaries(&nodes, &tets, cx);
            let v = initial_molar_volume();
            let moles = cavity.reference_volume_m3() / v;
            let mut gas = real_parcel(cavity.reference_volume_m3(), datum);
            close(gas.mass_kg(), MOLAR_MASS * moles, 2e-15);
            let heating = gas.add_heat(
                moles * energy(300.0, v, 360.0, v), solve_heat().gas, || false,
            ).unwrap();
            close(gas.state().temperature_k(), 360.0, 2e-7);
            close(heating.stored_energy_change_j, moles * energy(300.0, v, 360.0, v), 2e-9);
            let original = gas.clone();
            let system = ThermalPorousGasSystem::try_new(
                solid, &thermal, vec![entry(Arc::new(cavity), gas)],
                links(&thermal, 0.4), cx,
            ).unwrap();
            let observer = Times::default();
            let next = system.advance(
                &thermal, &observer, 1.0, &affine(&nodes, 1.0), &[], [0.0; 3], None,
                cavity_settings(), controls(), thermal_config(), config(),
                CouplingMethod::RelaxedPicard, &[0.02], &[0.0], solve_heat(), cx,
            ).unwrap();
            let decay: f64 = next.observer.0.iter().map(|dt| 1.0 - 0.4 * dt).product();
            let capacity: f64 = thermal.cells().values().map(|cell| cell.mass_kg * CP).sum();
            // Independent nonlinear finite-gas caloric balance; density-
            // dependent PR storage participates even at fixed gas volume.
            let (mut lo, mut hi) = (200.0, 600.0);
            for _ in 0..100 {
                let t = 0.5 * (lo + hi);
                let balance = moles * energy(360.0, v, t, v)
                    + capacity * (1.0 - decay) * (t - 300.0) - 0.02;
                if balance < 0.0 { lo = t; } else { hi = t; }
            }
            let expected_t = 0.5 * (lo + hi);
            let expected_wall = 300.0 + (1.0 - decay) * (expected_t - 300.0);
            let expected = pr(expected_t, v);
            let gas = &next.system.cavities()[0].gas;
            close(gas.state().temperature_k(), expected_t, 3e-7);
            close(gas.state().pressure_pa(), expected.pressure, 0.02);
            let bulk = gas.state().density_kg_m3() * gas.state().sound_speed_m_s().powi(2);
            close(bulk, expected.bulk, 0.05);
            assert!((bulk / (gas.state().gamma() * expected.pressure) - 1.0).abs() > 0.05);
            close(next.balance.total_gas_energy_change_j,
                moles * energy(360.0, v, expected_t, v), 2e-8);
            close(next.balance.gas_to_solid_heat_j,
                capacity * (expected_wall - 300.0), 2e-8);
            close(next.balance.external_solid_heat_j, 0.0, 1e-9);
            close(next.balance.combined_residual_j, 0.0, solve_heat().energy_tolerance_j);
            for cell in thermal.cells().values() {
                close(next.thermal.temperature_k(cell.id).unwrap(), expected_wall, 3e-7);
            }
            assert_eq!(gas.volume_m3().to_bits(), original.volume_m3().to_bits());
            assert_eq!(gas.mass_kg().to_bits(), original.mass_kg().to_bits());
            assert_eq!(gas.model(), original.model());
            assert_eq!(gas.point().axes().get("material-lot"), Some(&47.0));
            assert_eq!(system.cavities()[0].gas, original);
            assert!(observer.0.is_empty());
            let observed = [
                gas.state().temperature_k(), gas.state().pressure_pa(),
                next.balance.total_gas_energy_change_j, next.balance.gas_to_solid_heat_j,
            ];
            if let Some(reference) = without_datum {
                for i in 0..4 { close(observed[i], reference[i], if i == 1 { 2e-5 } else { 2e-10 }); }
            } else {
                without_datum = Some(observed);
            }
        }
    });
}

#[test]
fn heated_real_gas_retains_density_dependent_storage_across_compression_epochs() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let (nodes, tets) = shell();
        let mut thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let gas = real_parcel(cavity.reference_volume_m3(), 1e24);
        let initial = gas.clone();
        let moles = gas.mass_kg() / MOLAR_MASS;
        let mut system = ThermalPorousGasSystem::try_new(
            solid, &thermal, vec![entry(Arc::new(cavity), gas)],
            links(&thermal, 0.005), cx,
        ).unwrap();
        for (time, stretch, supplied) in [(1.0, 0.94, 0.015), (2.0, 0.90, -0.005)] {
            let before = system.cavities()[0].gas.clone();
            let old_cells = thermal.cells().clone();
            let old_states = system.solid().mechanics().element_states().to_vec();
            let supports = affine(&nodes, stretch);
            let refused = system.advance(
                &thermal, &(), time, &supports, &[], [0.0; 3], None,
                cavity_settings(), controls(), thermal_config(), config(),
                CouplingMethod::RelaxedPicard, &[supplied], &[1e12], solve_heat(), cx,
            );
            assert!(matches!(refused, Err(CavityHeatError::Gas(_))));
            assert_eq!(thermal.cells(), &old_cells);
            assert_eq!(system.cavities()[0].gas, before);
            assert_eq!(system.solid().mechanics().element_states(), old_states);
            let next = system.advance(
                &thermal, &(), time, &supports, &[], [0.0; 3], None,
                cavity_settings(), controls(), thermal_config(), config(),
                CouplingMethod::RelaxedPicard, &[supplied], &[0.0], solve_heat(), cx,
            ).unwrap();
            let after = &next.system.cavities()[0].gas;
            let t0 = before.state().temperature_k();
            let v0 = before.volume_m3() / moles;
            let t1 = after.state().temperature_k();
            let v1 = after.volume_m3() / moles;
            let stage_t = next.balance.gases[0].heat_stage_temperature_k;
            close(entropy(stage_t, v0, t1, v1), 0.0, 2e-8);
            close(after.state().pressure_pa(), pr(t1, v1).pressure, 0.05);
            close(next.balance.gases[0].rigid_heat_stage_energy_change_j,
                moles * energy(t0, v0, stage_t, v0), 3e-8);
            close(next.balance.gases[0].adiabatic_energy_change_j,
                moles * energy(stage_t, v0, t1, v1), 3e-8);
            close(next.balance.total_gas_energy_change_j,
                moles * energy(t0, v0, t1, v1), 5e-8);
            // Actual density departure participates in compression work.
            assert!(moles * (pr(t1, v1).residual_u - pr(stage_t, v0).residual_u).abs() > 1e-4);
            close(next.balance.gases[0].heat_residual_j, 0.0, solve_heat().energy_tolerance_j);
            close(next.balance.external_solid_heat_j, 0.0, 2e-8);
            close(next.balance.combined_residual_j, 0.0, solve_heat().energy_tolerance_j);
            assert_eq!(after.mass_kg().to_bits(), initial.mass_kg().to_bits());
            assert_eq!(after.model(), initial.model());
            let physical_volume = next.system.cavities()[0].cavity.evaluate(
                next.system.solid().mechanics().displacement_m().as_flattened(), cx,
            ).unwrap().volume_m3();
            assert_eq!(after.volume_m3().to_bits(), physical_volume.to_bits());
            for (old, new) in old_states.iter().zip(next.system.solid().mechanics().element_states()) {
                assert!(new.alpha >= old.alpha);
            }
            thermal = next.thermal;
            system = next.system;
        }
        assert_eq!(system.time_s(), 2.0);
        assert!(system.solid().mechanics().element_states().iter().all(|state| state.alpha > 0.0));
    });
}
