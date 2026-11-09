//! Native eddy viscosity drives original gas inventories and NASA-9 calorics.
//! The authored source fixture is algebraic, not measured-gas validation.

#[path = "nsf3_thermochemical/source.rs"]
mod source;

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    DensityProjectionControls, FlowFace, FlowField, FluidProperties, FvBoundary, FvFlow,
    FvSgsConfig, MaterialAdvectionControls, MaterialAdvectionInput, MaterialAdvectionStep,
    MaterialVolumeSource, SimpleConfig, UnsteadyConfig, VoxelDomain, simple_flow,
    simple_unsteady_with_sgs,
};
use fs_lbm::sgs3::{SgsModel, SgsUse};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::SpeciesNetworkSolve;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use std::collections::BTreeMap;

const R: f64 = 8.314_462_618_153_24;
const MOLAR_MASS: [f64; 2] = [0.028, 0.018_015_28];

fn energy() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 700.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}

fn controls() -> SpeciesNetworkSolve {
    SpeciesNetworkSolve {
        maximum_species_fraction: 0.2,
        mass_tolerance_kg: 1e-12,
        energy_tolerance_j: 1e-7,
    }
}

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1000.0,
        conductivity_w_m_k: 0.03,
        kinematic_viscosity_m2_s: 0.02,
    }
}

fn flow(domain: &VoxelDomain, eddy: &[f64]) -> FvFlow {
    let mut flow = simple_flow(
        domain,
        &fluid(),
        &SimpleConfig::new([FvBoundary::wall(); 6]),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    flow.eddy_viscosity_m2_s.copy_from_slice(eddy);
    flow
}

fn parcel(
    domain: &VoxelDomain,
    cell: usize,
    water: f64,
    t: f64,
    reference: f64,
) -> ClosedGasParcel {
    ClosedGasParcel::new(
        source::moist_model(water, reference),
        source::point(t, 100_000.0)
            .with("cell-lot", cell as f64)
            .unwrap(),
        domain.volume(cell),
    )
    .unwrap()
}

fn pair(domain: &VoxelDomain, reference: f64) -> BTreeMap<u64, ClosedGasParcel> {
    BTreeMap::from([
        (0, parcel(domain, 0, 0.1, 500.0, reference)),
        (1, parcel(domain, 1, 0.6, 350.0, reference)),
    ])
}

fn carrier(
    domain: &VoxelDomain,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
) -> MaterialAdvectionStep {
    let mass: Vec<_> = (0..domain.cell_count())
        .map(|cell| {
            parcels
                .get(&(cell as u64))
                .map_or(0.0, ClosedGasParcel::mass_kg)
        })
        .collect();
    // This accepted carrier binds mass and native geometry only. The original
    // parcels own all calorics in diffuse_turbulent_gas_species; a zero passive
    // advection scalar keeps its absolute-datum arithmetic out of that test,
    // including the deliberately enormous formation-reference shift below.
    let enthalpy = vec![0.0; domain.cell_count()];
    FlowField::quiescent(domain)
        .project_and_advect_material(
            domain,
            [FlowFace::Wall; 6],
            MaterialAdvectionInput {
                mass_kg: &mass,
                specific_enthalpy_j_kg: &enthalpy,
                volume_sources: &vec![MaterialVolumeSource::None; domain.cell_count()],
                boundary_inflow: [None; 6],
            },
            MaterialAdvectionControls {
                projection: DensityProjectionControls {
                    duration_s: 0.01,
                    relative_tolerance: 1e-12,
                    maximum_volume_defect_m3_s: 1e-10,
                    maximum_energy_defect_j: 1e-8,
                },
                max_cells: domain.cell_count(),
                max_connections: 6 * domain.cell_count(),
                max_iterations: 200,
                relative_tolerance: 1e-12,
                maximum_mass_defect_kg: 1e-12,
                maximum_enthalpy_defect_j: 1e-7,
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}

fn endpoints(parcels: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    parcels.keys().map(|&id| (id, energy())).collect()
}

fn fractions(parcel: &ClosedGasParcel) -> [f64; 2] {
    let x = parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    let denominator = x[0] * MOLAR_MASS[0] + x[1] * MOLAR_MASS[1];
    [
        x[0] * MOLAR_MASS[0] / denominator,
        x[1] * MOLAR_MASS[1] / denominator,
    ]
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn graded_turbulent_counterflow_carries_original_partial_enthalpy_and_reconstructs_temperature() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let parcels = pair(&domain, 0.0);
    let accepted = carrier(&domain, &parcels);
    let field = flow(&domain, &[0.0001, 0.0002]);
    let duration = 0.2;
    let schmidt = 0.7;
    let result = field
        .diffuse_turbulent_gas_species(
            &domain,
            &accepted,
            &parcels,
            schmidt,
            duration,
            &endpoints(&parcels),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.conducting_faces, 1);
    assert_eq!(result.exchange.step().fluxes.len(), 2);
    let y = [fractions(&parcels[&0]), fractions(&parcels[&1])];
    let mass = [parcels[&0].mass_kg(), parcels[&1].mass_kg()];
    let rho = [mass[0] / domain.volume(0), mass[1] / domain.volume(1)];
    let conductance =
        0.01 / (0.025 / (rho[0] * 0.0001 / schmidt) + 0.075 / (rho[1] * 0.0002 / schmidt));
    let transferred: [f64; 2] =
        core::array::from_fn(|k| duration * conductance * (y[0][k] - y[1][k]));
    let temperatures = [500.0, 350.0];
    let mut enthalpy = 0.0;
    for (k, name) in ["manufactured", "water"].into_iter().enumerate() {
        let row = result
            .exchange
            .step()
            .fluxes
            .iter()
            .find(|row| row.transfer.species.as_str() == name)
            .unwrap();
        let signed = if row.transfer.donor == 0 {
            row.transfer.mass_kg
        } else {
            -row.transfer.mass_kg
        };
        near(signed, transferred[k], 2e-17);
        let donor_t = temperatures[usize::from(transferred[k] < 0.0)];
        // Independent integral of the authored cp/R = 3 + 0.001 T.
        let h = R / MOLAR_MASS[k] * (3.0 * donor_t + 0.0005 * donor_t * donor_t);
        near(row.source_enthalpy_j_kg, h, 1e-8);
        enthalpy += transferred[k] * h;
    }
    for (cell, sign) in [(0, -1.0), (1, 1.0)] {
        let next = &result.exchange.parcels()[&(cell as u64)];
        let initial_mass = y[cell].map(|fraction| fraction * mass[cell]);
        let final_mass: [f64; 2] =
            core::array::from_fn(|k| initial_mass[k] + sign * transferred[k]);
        let capacity = |m: [f64; 2]| (0..2).map(|k| m[k] * R / MOLAR_MASS[k]).sum::<f64>();
        let t = temperatures[cell];
        let initial_u = capacity(initial_mass) * (2.0 * t + 0.0005 * t * t);
        let target = (initial_u + sign * enthalpy) / capacity(final_mass);
        let expected_t = 2.0 * target / (2.0 + (4.0 + 0.002 * target).sqrt());
        near(next.state().state().temperature, expected_t, 2e-7);
        near(next.mass_kg(), mass[cell], controls().mass_tolerance_kg);
        let new_y = fractions(next);
        for k in 0..2 {
            near(new_y[k] * next.mass_kg(), final_mass[k], 2e-16);
        }
        near(
            next.state().state().pressure,
            next.mass_kg() / next.volume_m3()
                * next.state().state().specific_gas_constant
                * expected_t,
            1e-4,
        );
        assert_eq!(next.point().axes().get("cell-lot"), Some(&(cell as f64)));
        assert_eq!(next.volume_m3(), domain.volume(cell));
        assert_eq!(
            next.model().source_packs().len(),
            parcels[&(cell as u64)].model().source_packs().len()
        );
        assert_ne!(next.state().state().temperature, t);
    }
    assert!(result.exchange.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(result.maximum_face_mass_defect_kg <= controls().mass_tolerance_kg);
    assert!(result.maximum_carrier_mass_defect_kg <= controls().mass_tolerance_kg);
}

#[test]
fn formation_reference_shift_leaves_the_material_endpoint_and_transport_unchanged() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let field = flow(&domain, &[0.0001, 0.0002]);
    let solve = |reference| {
        let parcels = pair(&domain, reference);
        field
            .diffuse_turbulent_gas_species(
                &domain,
                &carrier(&domain, &parcels),
                &parcels,
                0.7,
                0.2,
                &endpoints(&parcels),
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let original = solve(0.0);
    let shifted = solve(1e18);
    for id in [0, 1] {
        let a = &original.exchange.parcels()[&id];
        let b = &shifted.exchange.parcels()[&id];
        near(
            a.state().state().temperature,
            b.state().state().temperature,
            1e-8,
        );
        near(a.state().state().pressure, b.state().state().pressure, 1e-5);
        near(a.mass_kg(), b.mass_kg(), 1e-16);
        assert_eq!(
            a.state().state().dynamic_viscosity,
            b.state().state().dynamic_viscosity
        );
        assert_eq!(
            a.state().state().thermal_conductivity,
            b.state().state().thermal_conductivity
        );
    }
    near(
        original.exchange.step().total_energy_residual_j,
        shifted.exchange.step().total_energy_residual_j,
        1e-10,
    );
}

#[test]
fn native_sgs_field_drives_source_owned_composition_and_caloric_updates() {
    let domain = VoxelDomain::new(3, 3, 2, 0.25).unwrap();
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall {
        velocity: [0.2, 0.0, 0.0],
    };
    let mut march = UnsteadyConfig::new(0.02, 2);
    march.inner_iterations = 400;
    march.inner_tolerance = 1e-8;
    let gate = CancelGate::new_clock_free();
    let field = simple_unsteady_with_sgs(
        &domain,
        &fluid(),
        &SimpleConfig::new(faces),
        &march,
        FvSgsConfig {
            model: SgsModel::Smagorinsky { cs: 0.2 },
            filter_width_factor: 1.0,
            maximum_eddy_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        |_| 1.0,
        &gate,
    )
    .unwrap()
    .flow;
    assert!(field.eddy_viscosity_m2_s.iter().any(|&nu| nu > 0.0));
    let parcels: BTreeMap<_, _> = (0..domain.cell_count())
        .map(|cell| {
            let x = domain.coords(cell)[0];
            (
                cell as u64,
                parcel(
                    &domain,
                    cell,
                    if x == 0 { 0.1 } else { 0.6 },
                    if x == 0 { 450.0 } else { 350.0 },
                    -4000.0,
                ),
            )
        })
        .collect();
    let accepted = carrier(&domain, &parcels);
    let run = |flow: &FvFlow, sc| {
        flow.diffuse_turbulent_gas_species(
            &domain,
            &accepted,
            &parcels,
            sc,
            0.01,
            &endpoints(&parcels),
            controls(),
            &gate,
        )
        .unwrap()
    };
    let result = run(&field, 0.7);
    assert!(!result.exchange.step().fluxes.is_empty());
    assert!(result.exchange.parcels().iter().any(|(id, p)| {
        p.model().thermodynamics().mole_composition()
            != parcels[id].model().thermodynamics().mole_composition()
    }));
    assert!(result.exchange.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    let mut scaled = field.clone();
    scaled
        .eddy_viscosity_m2_s
        .iter_mut()
        .for_each(|nu| *nu *= 2.0);
    scaled.heat_eddy_viscosity_m2_s.fill(100.0);
    scaled.turbulent_prandtl = 9.0;
    assert_eq!(result, run(&scaled, 1.4));
}

#[test]
fn absent_species_enter_through_original_sources_and_solid_faces_stay_closed() {
    use fs_lbm::conjugate::Voxel;
    use fs_qty::Amount;

    let mut domain = VoxelDomain::new(3, 1, 1, 0.1).unwrap();
    domain.set(2, 0, 0, Voxel::Solid(0));
    let catalog = source::moist_model(0.5, -4000.0);
    let mut parcels = BTreeMap::new();
    for cell in 0..2 {
        let amounts = if cell == 0 { [1.0, 0.0] } else { [0.0, 1.0] };
        let model = catalog
            .with_mole_amounts(&amounts.map(Amount::new))
            .unwrap();
        assert_eq!(model.thermodynamics().components().len(), 1);
        parcels.insert(
            cell as u64,
            ClosedGasParcel::new(model, source::point(400.0, 100_000.0), domain.volume(cell))
                .unwrap(),
        );
    }
    let accepted = carrier(&domain, &parcels);
    let field = flow(&domain, &[0.001, 0.001, 0.0]);
    let result = field
        .diffuse_turbulent_gas_species(
            &domain,
            &accepted,
            &parcels,
            0.7,
            0.05,
            &endpoints(&parcels),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.conducting_faces, 1);
    assert_eq!(result.exchange.step().fluxes.len(), 2);
    assert_eq!(result.exchange.parcels().len(), 2);
    assert!(!result.exchange.parcels().contains_key(&2));
    for id in [0, 1] {
        let original = &parcels[&id];
        let next = &result.exchange.parcels()[&id];
        assert_eq!(next.model().thermodynamics().components().len(), 2);
        let incoming = result
            .exchange
            .step()
            .fluxes
            .iter()
            .find(|row| row.transfer.receiver == id)
            .unwrap();
        let species = next
            .model()
            .thermodynamics()
            .components()
            .iter()
            .position(|s| s.species() == &incoming.transfer.species)
            .unwrap();
        near(
            fractions(next)[species] * next.mass_kg(),
            incoming.transfer.mass_kg,
            1e-16,
        );
        near(
            next.mass_kg(),
            original.mass_kg(),
            controls().mass_tolerance_kg,
        );
        assert_eq!(original.model().thermodynamics().components().len(), 1);
    }
    assert!(result.exchange.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
}

#[test]
fn disabled_and_uniform_composition_are_exact_and_refusal_changes_no_owner() {
    let domain = VoxelDomain::new(2, 1, 1, 0.1).unwrap();
    let parcels = pair(&domain, 0.0);
    let accepted = carrier(&domain, &parcels);
    let saved = (parcels.clone(), accepted.clone());
    let field = flow(&domain, &[0.0; 2]);
    let gate = CancelGate::new_clock_free();
    let disabled = field
        .diffuse_turbulent_gas_species(
            &domain,
            &accepted,
            &parcels,
            0.7,
            10.0,
            &endpoints(&parcels),
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(disabled.exchange.parcels(), &parcels);
    assert!(disabled.exchange.step().fluxes.is_empty());
    let uniform = BTreeMap::from([
        (0, parcel(&domain, 0, 0.25, 500.0, 0.0)),
        (1, parcel(&domain, 1, 0.25, 350.0, 0.0)),
    ]);
    let turbulent = flow(&domain, &[0.01; 2]);
    let unchanged = turbulent
        .diffuse_turbulent_gas_species(
            &domain,
            &carrier(&domain, &uniform),
            &uniform,
            0.7,
            10.0,
            &endpoints(&uniform),
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(unchanged.exchange.parcels(), &uniform);
    assert!(unchanged.exchange.step().fluxes.is_empty());
    for sc in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            turbulent
                .diffuse_turbulent_gas_species(
                    &domain,
                    &accepted,
                    &parcels,
                    sc,
                    0.1,
                    &endpoints(&parcels),
                    controls(),
                    &gate,
                )
                .is_err()
        );
    }
    assert!(
        turbulent
            .diffuse_turbulent_gas_species(
                &domain,
                &accepted,
                &parcels,
                0.7,
                1e6,
                &endpoints(&parcels),
                controls(),
                &gate,
            )
            .is_err()
    );
    let mut changed = accepted.clone();
    changed.mass_kg[0] *= 1.01;
    assert!(
        turbulent
            .diffuse_turbulent_gas_species(
                &domain,
                &changed,
                &parcels,
                0.7,
                0.1,
                &endpoints(&parcels),
                controls(),
                &gate,
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        turbulent
            .diffuse_turbulent_gas_species(
                &domain,
                &accepted,
                &parcels,
                0.7,
                0.1,
                &endpoints(&parcels),
                controls(),
                &cancelled,
            )
            .is_err()
    );
    assert_eq!((parcels, accepted), saved);
}
