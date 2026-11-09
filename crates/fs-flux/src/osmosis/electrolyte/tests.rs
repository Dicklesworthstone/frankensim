use super::*;
use fs_material::osmosis::{
    HydraulicPermeance, ReflectionCoefficient, SolutePermeance,
    electrodiffusion::ElectrodiffusiveMembrane, solubility::electrolyte::SolventActivityBasis,
};
use fs_qty::{AmountConcentration, Volume};
use fs_thermochem::{
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    stoichiometry::thermodynamics::acid_base::{
        SpeciationOptions, SpectatorIon,
        precipitation::kinetics::{MineralGrowthPolicyV1, MineralKineticsV1, MineralSurfaceRateV1},
    },
};

// Reuse the original source-qualified authored chemistry, not a replacement
// speciation or source model. Its numerical values are test data, not measured.
#[path = "../../../../fs-material/src/osmosis/solubility/electrolyte/tests/fixture.rs"]
mod fixture;
use fixture::*;

const T: f64 = 298.15;
const V: f64 = 1e-3;
const RATE: f64 = 1e-3;

fn aqueous_options() -> SpeciationOptions {
    SpeciationOptions {
        absolute_charge_tolerance: Amount::new(1e-18),
        relative_charge_tolerance: Dimensionless::new(1e-12),
        ph_tolerance: Dimensionless::new(1e-13),
        max_iterations: 256,
        ..SpeciationOptions::default()
    }
}
fn phase_controls() -> MineralKineticControlsV1 {
    MineralKineticControlsV1 {
        aqueous: aqueous_options(),
        initial_step: Time::new(0.02),
        minimum_step: Time::new(1e-9),
        maximum_step: Time::new(0.1),
        relative_tolerance: Dimensionless::new(1e-7),
        amount_absolute_tolerance: Amount::new(1e-12),
        extent_absolute_tolerance: Amount::new(1e-13),
        maximum_stage_iterations: 128,
        maximum_attempts: 8192,
    }
}
fn outer(relative: f64, step: f64) -> OsmoticStepControls {
    OsmoticStepControls {
        maximum_step_s: step,
        relative_tolerance: relative,
        absolute_volume_tolerance_m3: 1e-12,
        absolute_amount_tolerance_mol: 1e-10,
        maximum_steps: 20_000,
        maximum_rejections: 128,
        ..OsmoticStepControls::default()
    }
}
fn salt_model(rate: f64) -> MineralKineticsV1 {
    MineralKineticsV1::from_salt(
        salt(T, 4e-6),
        MineralSurfaceRateV1::new(rate, MineralGrowthPolicyV1::ExistingSolidOnly).unwrap(),
    )
    .unwrap()
}
fn material(model: MineralKineticsV1, dissolved: f64) -> KineticElectrolyteSolubilitySolution {
    KineticElectrolyteSolubilitySolution::new(
        model,
        SolventActivityBasis::new(18e-6, AmountConcentration::new(1000.0)).unwrap(),
        Volume::new(V),
        Amount::new(dissolved),
        Amount::new(0.11 - dissolved),
        vec![(FAMILY.into(), Amount::new(0.02 + 2.0 * dissolved))],
        vec![],
        aqueous_options(),
        || false,
    )
    .unwrap()
}
fn cell(material: KineticElectrolyteSolubilitySolution, pressure: f64) -> OsmoticCompartment {
    OsmoticCompartment::with_kinetic_electrolyte_solubility(
        material,
        phase_controls(),
        OsmoticPressureVolume::new(V, pressure, 0.0).unwrap(),
        1e-5,
    )
    .unwrap()
}
fn membrane(hydraulic: f64, reflection: f64, permeance: f64, count: usize) -> OsmoticMembrane {
    OsmoticMembrane::new(
        0,
        1,
        1.0,
        MembraneTransport::new(
            HydraulicPermeance::new(hydraulic).unwrap(),
            vec![
                (
                    ReflectionCoefficient::new(reflection).unwrap(),
                    SolutePermeance::new(permeance).unwrap()
                );
                count
            ],
        )
        .unwrap(),
    )
    .unwrap()
}
fn cells(dissolved: f64) -> Vec<OsmoticCompartment> {
    let model = salt_model(RATE);
    vec![
        cell(material(model.clone(), dissolved), 50_000.0),
        cell(material(model, dissolved), 0.0),
    ]
}
fn axis(cells: &[OsmoticCompartment]) -> Vec<String> {
    cells[0]
        .kinetic_electrolyte
        .as_ref()
        .unwrap()
        .material
        .species()
        .iter()
        .map(|id| id.as_str().into())
        .collect()
}
fn network(dissolved: f64, hydraulic: f64) -> OsmoticNetwork {
    let cells = cells(dissolved);
    let names = axis(&cells);
    let membrane = membrane(hydraulic, 1.0, 0.0, names.len());
    OsmoticNetwork::new(names, cells, vec![membrane], 1 << 22).unwrap()
}
fn phase(cell: &OsmoticCompartment) -> &KineticElectrolyteSolubilitySolution {
    &cell.kinetic_electrolyte.as_ref().unwrap().material
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "actual={actual:.16e} expected={expected:.16e} tolerance={absolute:e}"
    );
}

// Independent continuous limit: solve the authored charge cubic in log(h),
// then integrate only V_left and the two dissolved formula inventories.
// No production aqueous/kinetic/water stage is called by this RK4 reference.
fn chemical(volume: f64, dissolved: f64) -> (f64, f64) {
    let reference = 1000.0 * volume;
    let family = 0.02 + 2.0 * dissolved;
    let (mut lo, mut hi) = (1e-14_f64.ln(), 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (lo + hi);
        let h = middle.exp();
        let charge = reference * (h - 1e-14 / h) + 2.0 * dissolved - family * 1e-7 / (h + 1e-7);
        if charge > 0.0 {
            hi = middle;
        } else {
            lo = middle;
        }
    }
    let h = (0.5 * (lo + hi)).exp();
    let anion = family * 1e-7 / (h + 1e-7);
    let omega = (dissolved / reference) * (anion / reference).powi(2) / 4e-6;
    let particles = dissolved + family + reference * (h + 1e-14 / h);
    (particles / volume, omega)
}
fn derivative(state: [f64; 3], hydraulic: f64) -> [f64; 3] {
    let (left, omega_l) = chemical(state[0], state[1]);
    let (right, omega_r) = chemical(2.0 * V - state[0], state[2]);
    [
        hydraulic * (R * T * (left - right) - 50_000.0),
        RATE * (1.0 - omega_l),
        RATE * (1.0 - omega_r),
    ]
}
fn reference(duration: f64, dt: f64, hydraulic: f64) -> [f64; 3] {
    let mut state = [V, 0.01, 0.01];
    let count = (duration / dt).round() as usize;
    let dt = duration / count as f64;
    for _ in 0..count {
        let a = derivative(state, hydraulic);
        let b = derivative(
            core::array::from_fn(|i| state[i] + 0.5 * dt * a[i]),
            hydraulic,
        );
        let c = derivative(
            core::array::from_fn(|i| state[i] + 0.5 * dt * b[i]),
            hydraulic,
        );
        let d = derivative(core::array::from_fn(|i| state[i] + dt * c[i]), hydraulic);
        state =
            core::array::from_fn(|i| state[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    state
}
fn observed(network: &OsmoticNetwork) -> [f64; 3] {
    [
        network.cells[0].solution.solvent_volume_m3(),
        phase(&network.cells[0]).dissolved_cation().value(),
        phase(&network.cells[1]).dissolved_cation().value(),
    ]
}

#[test]
fn water_dilution_and_concentration_drive_original_salt_kinetics_and_refine() {
    let hydraulic = 1e-9;
    let expected = reference(1.0, 0.0005, hydraulic);
    let finer_reference = reference(1.0, 0.00025, hydraulic);
    for (a, b) in expected.into_iter().zip(finer_reference) {
        near(a, b, 2e-13);
    }
    let mut coarse = network(0.01, hydraulic);
    let mut fine = coarse.clone();
    let old_source = phase(&fine.cells[0]).aqueous().model().water() as *const _;
    coarse.advance(1.0, outer(1e-3, 0.2), || false).unwrap();
    // The water stage is backward Euler: local step-doubling tolerance is
    // not a global error bound. Refine dt for the independent finite-window
    // physical oracle while retaining the same local/native accuracy.
    let report = fine.advance(1.0, outer(1e-5, 0.005), || false).unwrap();
    let actual = observed(&fine);
    assert!(actual[0] < V && fine.cells[1].solution.solvent_volume_m3() > V);
    assert!(actual[1] < 0.01 && actual[2] > 0.01);
    assert!(phase(&fine.cells[0]).solid_amount().value() > 0.1);
    assert!(phase(&fine.cells[1]).solid_amount().value() < 0.1);
    let mut fine_error = 0.0;
    let mut coarse_error = 0.0;
    for i in 0..3 {
        near(actual[i], expected[i], if i == 0 { 3e-8 } else { 3e-7 });
        fine_error += ((actual[i] - expected[i]) / expected[i]).abs();
        coarse_error += ((observed(&coarse)[i] - expected[i]) / expected[i]).abs();
    }
    assert!(
        fine_error < coarse_error,
        "fine={fine_error:e}, coarse={coarse_error:e}"
    );
    assert_eq!(
        old_source,
        phase(&fine.cells[0]).aqueous().model().water() as *const _
    );
    assert!(report.maximum_scaled_error <= 1.0);
    assert!(report.phase_evaluations > 0 && report.electrolyte_volume_speciations > 0);
    near(report.water_balance_residual_m3, 0.0, 2e-18);
    assert_eq!(report.species_count, 0);
    assert_eq!(
        report.solute_balance_residual_mol,
        [0.0; MAX_OSMOTIC_SPECIES]
    );
    let balances = report.electrolyte_balances.unwrap();
    assert_eq!(balances.len(), 2);
    for balance in balances {
        assert_eq!(balance.inventories.len(), 2);
        assert!(matches!(
            balance.inventories[0].inventory,
            ElectrolyteInventory::Cation(_)
        ));
        assert!(matches!(
            balance.inventories[1].inventory,
            ElectrolyteInventory::Family(_)
        ));
        for inventory in balance.inventories {
            // Match the original mineral owner's representability allowance
            // for actual independent phase stocks. This is floating roundoff
            // after many accepted stages, not a temporal-error allowance.
            near(
                inventory.residual_mol,
                0.0,
                64.0 * f64::EPSILON * inventory.total_mol,
            );
        }
        near(balance.charge_balance_residual_mol, 0.0, 1e-12);
    }
    assert_eq!(
        phase(&fine.cells[0]).component_total(),
        phase(&coarse.cells[0]).component_total()
    );
    fine.advance(0.1, outer(1e-5, 0.05), || false).unwrap();
    assert!(phase(&fine.cells[0]).solid_amount().value() > 0.1);
}

#[test]
fn sealed_water_recovers_original_finite_rate_mineral_evolution() {
    let mut network = network(0.003, 0.0);
    let original = phase(&network.cells[0]).clone();
    let direct = original
        .advance(Time::new(0.2), phase_controls(), || false)
        .unwrap();
    let report = network.advance(0.2, outer(1e-5, 0.05), || false).unwrap();
    assert_eq!(network.cells[0].solution.solvent_volume_m3(), V);
    assert_eq!(network.cells[1].solution.solvent_volume_m3(), V);
    near(
        phase(&network.cells[0]).dissolved_cation().value(),
        direct.solution().dissolved_cation().value(),
        3e-9,
    );
    near(
        phase(&network.cells[0]).solid_amount().value(),
        direct.solution().solid_amount().value(),
        3e-9,
    );
    assert_eq!(report.electrolyte_volume_speciations, 0);
    assert!(report.phase_evaluations > 0);
}

#[test]
fn default_outer_accuracy_retains_a_resolvable_relative_phase_allowance() {
    let mut coupled = network(0.01, 1e-9);
    let report = coupled
        .advance(0.1, OsmoticStepControls::default(), || false)
        .unwrap();
    assert!(coupled.cells[0].solution.solvent_volume_m3() < V);
    assert!(phase(&coupled.cells[0]).solid_amount().value() > 0.1);
    assert!(phase(&coupled.cells[1]).solid_amount().value() < 0.1);
    assert!(report.maximum_scaled_error <= 1.0);
    assert!(report.electrolyte_balances.is_some());

    // A fresh solid-bearing bath must also admit the original exact-zero
    // dissolved boundary; no dissolved concentration floor is inserted.
    let mut fresh = network(0.0, 0.0);
    let report = fresh
        .advance(0.1, OsmoticStepControls::default(), || false)
        .unwrap();
    assert!(phase(&fresh.cells[0]).dissolved_cation().value() > 0.0);
    assert!(phase(&fresh.cells[0]).solid_amount().value() < 0.11);
    assert_eq!(fresh.cells[0].solution.solvent_volume_m3(), V);
    assert!(report.maximum_scaled_error <= 1.0);

    // The explicitly selected seeded-surface law can create its first actual
    // solid. The new phase's relative scale belongs in the root allowance.
    let source = MineralKineticsV1::from_salt(
        salt(T, 4e-6),
        MineralSurfaceRateV1::new(RATE, MineralGrowthPolicyV1::SeededSurface).unwrap(),
    )
    .unwrap();
    let make = |model| {
        KineticElectrolyteSolubilitySolution::new(
            model,
            SolventActivityBasis::new(18e-6, AmountConcentration::new(1000.0)).unwrap(),
            Volume::new(V),
            Amount::new(0.02),
            Amount::new(0.0),
            vec![(FAMILY.into(), Amount::new(0.06))],
            vec![],
            aqueous_options(),
            || false,
        )
        .unwrap()
    };
    let seeded_cells = vec![cell(make(source.clone()), 0.0), cell(make(source), 0.0)];
    let names = axis(&seeded_cells);
    let face = membrane(0.0, 1.0, 0.0, names.len());
    let mut seeded = OsmoticNetwork::new(names, seeded_cells, vec![face], 1 << 22).unwrap();
    let report = seeded
        .advance(0.02, OsmoticStepControls::default(), || false)
        .unwrap();
    assert!(phase(&seeded.cells[0]).solid_amount().value() > 0.0);
    assert!(phase(&seeded.cells[0]).dissolved_cation().value() < 0.02);
    assert_eq!(phase(&seeded.cells[0]).component_total(), Amount::new(0.02));
    assert!(report.maximum_scaled_error <= 1.0);
}

#[test]
fn hydrolysis_water_candidates_freeze_mineral_and_resolve_p_h() {
    let model = MineralKineticsV1::from_hydrolysis(
        hydrolysis(T, 1e-9),
        MineralSurfaceRateV1::new(1e-6, MineralGrowthPolicyV1::ExistingSolidOnly).unwrap(),
    )
    .unwrap();
    let make = |model| {
        KineticElectrolyteSolubilitySolution::new(
            model,
            SolventActivityBasis::new(18e-6, AmountConcentration::new(1000.0)).unwrap(),
            Volume::new(V),
            Amount::new(1e-5),
            Amount::new(1e-4),
            vec![],
            vec![SpectatorIon {
                species: id("Cl[-](aq)"),
                charge: -1,
                amount: Amount::new(3e-5),
            }],
            aqueous_options(),
            || false,
        )
        .unwrap()
    };
    let left = cell(make(model.clone()), 0.0);
    let right = cell(make(model), 0.0);
    let face = membrane(1e-9, 1.0, 0.0, left.solution.solute_moles().len());
    let mut calls = 0;
    let endpoint =
        pair_endpoint(&left, &right, &face, 1.0, 2e-4, &mut calls, &mut || false).unwrap();
    let [a, b] = endpoint.electrolytes.unwrap();
    assert_eq!(calls, 2);
    for material in [&a, &b] {
        assert_eq!(material.solid_amount(), Amount::new(1e-4));
        assert_eq!(material.dissolved_cation(), Amount::new(1e-5));
        assert_eq!(material.component_total(), phase(&left).component_total());
    }
    assert!(a.aqueous().ph.value() < phase(&left).aqueous().ph.value() - 0.09);
    assert!(b.aqueous().ph.value() > phase(&right).aqueous().ph.value() + 0.07);
    let cells = vec![left, right];
    let mut network = OsmoticNetwork::new(axis(&cells), cells, vec![face], 1 << 22).unwrap();
    let report = network.advance(0.05, outer(1e-4, 0.02), || false).unwrap();
    assert!(phase(&network.cells[0]).dissolved_cation().value() > 1e-5);
    assert!(
        report.electrolyte_balances.unwrap()[0]
            .inventories
            .iter()
            .any(|entry| matches!(entry.inventory, ElectrolyteInventory::Spectator { .. }))
    );
}

fn snapshot(network: &OsmoticNetwork) -> Vec<u64> {
    let mut values = vec![network.time_s().to_bits()];
    for cell in &network.cells {
        let material = phase(cell);
        values.extend([
            cell.solution.solvent_volume_m3().to_bits(),
            material.solid_amount().value().to_bits(),
            material.dissolved_cation().value().to_bits(),
            material.aqueous().ph.value().to_bits(),
        ]);
        values.extend(cell.solution.solute_moles().iter().map(|n| n.to_bits()));
    }
    values
}

#[test]
fn cancellation_and_work_refusal_publish_no_ionic_or_water_state() {
    let original = network(0.01, 1e-9);
    let mut completed = original.clone();
    let controls = outer(1e-4, 0.02);
    let mut polls = 0usize;
    completed
        .advance(0.03, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut cancelled = original.clone();
    let before = snapshot(&cancelled);
    let mut seen = 0usize;
    assert!(matches!(
        cancelled.advance(0.03, controls, || {
            seen += 1;
            seen >= polls - 1
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(snapshot(&cancelled), before);
    let mut limited = original;
    let mut controls = controls;
    controls.maximum_steps = 1;
    controls.maximum_step_s = 0.005;
    assert!(matches!(
        limited.advance(0.03, controls, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(snapshot(&limited), before);
}

#[test]
fn unsupported_ionic_permeation_diffusion_and_foreign_sources_refuse() {
    let cells = cells(0.01);
    let names = axis(&cells);
    for (reflection, permeance) in [(1.0, 1e-7), (0.99, 0.0)] {
        assert!(
            OsmoticNetwork::new(
                names.clone(),
                cells.clone(),
                vec![membrane(1e-9, reflection, permeance, names.len())],
                1 << 22
            )
            .is_err()
        );
    }
    let mut foreign = cells.clone();
    foreign[1] = cell(material(salt_model(RATE), 0.01), 0.0);
    assert!(OsmoticNetwork::new(names.clone(), foreign, vec![], 1 << 22).is_err());
    assert!(
        OsmoticNetwork::with_solution_diffusion(
            names.clone(),
            cells.clone(),
            vec![],
            diffusion::SolutionDiffusionControls {
                constitutive: fs_material::osmosis::diffusion::SolutionDiffusionSolve {
                    friction: fs_thermochem::transport::maxwell_stefan::MaxwellStefanMolarSolveV1::default(),
                    entropy_tolerance_w_m3_k: 1e-9,
                },
                maximum_species_fraction: 0.2,
                amount_balance_tolerance_mol: 1e-16,
                free_energy_tolerance_j: 1e-12,
            },
            1 << 22
        )
        .is_err()
    );
    let electrical = ElectrodiffusiveMembrane::new(
        HydraulicPermeance::new(1e-9).unwrap(),
        vec![1, -1, 0, -1, 2],
        vec![SolutePermeance::new(0.0).unwrap(); names.len()],
    )
    .unwrap();
    assert!(
        electrodiffusion::IonicOsmoticPair::new(
            cells.try_into().unwrap(),
            electrical,
            [0.0; 2],
            [0.0; 2],
            1.0,
            1e-12,
            1 << 22
        )
        .is_err()
    );
}
