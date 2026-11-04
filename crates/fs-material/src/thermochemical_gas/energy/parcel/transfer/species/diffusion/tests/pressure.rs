//! Independent binary barodiffusion and conservative network integration.
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::diffusion::network as molecular;
use crate::thermochemical_gas::energy::parcel::transfer::species::network::SpeciesNetworkSolve;
use fs_qty::Pressure;

fn driving() -> GasDiffusionDriving {
    GasDiffusionDriving::CompositionAndPressureMass { mass_flux_tolerance_kg_m2_s: 1e-16 }
}
fn pressure_controls() -> GasDiffusionSolve {
    GasDiffusionSolve { maximum_pressure_difference_pa: 2e5,
        maximum_temperature_difference_k: 1.0, ..controls() }
}
fn gas_at(x: f64, pressure: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], 0.001), point(400.0, pressure), 0.001).unwrap()
}

#[test]
fn pressure_driving_matches_binary_oracle_and_full_chemical_entropy() {
    let (a, b) = (gas_at(0.5, 120000.0), gas_at(0.5, 80000.0));
    let d = card(2e-5); let context = point(300.0, 90000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let old = a.diffusion_flux_to(&b, geometry(), &pairs, pressure_controls(), || false).unwrap();
    assert!(old.species.iter().all(|row| row.molar_flux_mol_m2_s == 0.0));
    let flux = a.diffusion_flux_to_with_driving(&b, geometry(), &pairs,
        pressure_controls(), driving(), || false).unwrap();
    assert_eq!(flux.driving, driving());
    // Independently eliminate J_B = -M_A/M_B J_A in the original binary row.
    let (ma, mb, x) = (0.028, 0.004, 0.5);
    let mean_mass = x*ma + (1.0-x)*mb;
    let gradient = (80000.0_f64/120000.0).ln()/geometry().distance_m;
    let force = (x-x*ma/mean_mass)*gradient;
    let c = 100000.0/(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*400.0);
    let ja = -c*2e-5*force*mb/mean_mass;
    let jb = c*2e-5*force*ma/mean_mass;
    near(flux.species[0].molar_flux_mol_m2_s, ja, 1e-14);
    near(flux.species[1].molar_flux_mol_m2_s, jb, 1e-14);
    near(flux.log_pressure_gradient_per_m, gradient, 1e-12);
    near(flux.mass_frame_residual_kg_m2_s, 0.0, 1e-18);
    near(flux.molar_frame_residual_mol_m2_s, ja+jb, 1e-14);
    assert!(flux.molar_frame_residual_mol_m2_s.abs() > 1e-3);
    // Uniform composition: full chemical force is grad ln(p), NOT zero.
    let entropy = -UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*(ja+jb)*gradient;
    assert!(entropy > 0.0);
    near(flux.drag_entropy_production_w_m3_k.unwrap(), entropy, 1e-11);
    let reverse = b.diffusion_flux_to_with_driving(&a, geometry(), &pairs,
        pressure_controls(), driving(), || false).unwrap();
    assert_eq!(reverse.log_pressure_gradient_per_m, -flux.log_pressure_gradient_per_m);
    for (left, right) in flux.species.iter().zip(&reverse.species) {
        near(left.molar_flux_mol_m2_s, -right.molar_flux_mol_m2_s, 1e-15);
        near(left.mass_rate_left_to_right_kg_s, -right.mass_rate_left_to_right_kg_s, 1e-18);
    }
    near(reverse.drag_entropy_production_w_m3_k.unwrap(), entropy, 1e-11);
}

#[test]
fn equal_pressure_is_the_correct_mass_frame_transform_of_molar_counterdiffusion() {
    let (a, b) = (gas(0.75), gas(0.25)); let d = card(2e-5);
    let context = point(400.0, 100000.0); let pairs = [pair("A", "B", &d, &context)];
    let molar = a.diffusion_flux_to(&b, geometry(), &pairs, controls(), || false).unwrap();
    let mass = a.diffusion_flux_to_with_driving(&b, geometry(), &pairs,
        controls(), driving(), || false).unwrap();
    assert_eq!(mass.log_pressure_gradient_per_m, 0.0);
    let correction = (0.028*molar.species[0].molar_flux_mol_m2_s
        + 0.004*molar.species[1].molar_flux_mol_m2_s)/0.016;
    for (a, b) in molar.species.iter().zip(&mass.species) {
        near(b.molar_flux_mol_m2_s, a.molar_flux_mol_m2_s-0.5*correction, 1e-14);
    }
    near(mass.mass_frame_residual_kg_m2_s, 0.0, 1e-18);
    near(mass.drag_entropy_production_w_m3_k.unwrap(),
        molar.drag_entropy_production_w_m3_k.unwrap(), 1e-11);
}

#[test]
fn equal_molar_masses_remove_pressure_separation_and_frame_difference() {
    let gas = |x, p| {
        let mixture = FrozenIdealGasMixtureModelV1::new(
            vec![species("A", 0.028, 3.5, 0.001), species("B", 0.028, 2.5, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0-x]).unwrap()).unwrap();
        let material = ThermochemicalGasModel::new(mixture,
            vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)], None).unwrap();
        ClosedGasParcel::new(material, point(400.0, p), 0.001).unwrap()
    };
    let d = card(2e-5); let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    for (xl, xr) in [(0.75, 0.25), (0.5, 0.5)] {
        let (a, b) = (gas(xl, 120000.0), gas(xr, 80000.0));
        let old = a.diffusion_flux_to(&b, geometry(), &pairs, pressure_controls(), || false).unwrap();
        let new = a.diffusion_flux_to_with_driving(&b, geometry(), &pairs,
            pressure_controls(), driving(), || false).unwrap();
        for (a, b) in old.species.iter().zip(&new.species) {
            near(a.molar_flux_mol_m2_s, b.molar_flux_mol_m2_s, 1e-14);
        }
    }
}

#[test]
fn adjacent_pressures_keep_the_small_driving_force() {
    let low = 100000.0_f64; let high = f64::from_bits(low.to_bits()+1);
    let (a, b) = (gas_at(0.5, low), gas_at(0.5, high));
    let d = card(2e-5); let context = point(400.0, low);
    let pairs = [pair("A", "B", &d, &context)];
    let flux = a.diffusion_flux_to_with_driving(&b, geometry(), &pairs,
        controls(), driving(), || false).unwrap();
    near(flux.log_pressure_gradient_per_m,
        ((high-low)/low).ln_1p()/geometry().distance_m, 1e-28);
    assert!(flux.log_pressure_gradient_per_m > 0.0);
    assert!(flux.species[0].molar_flux_mol_m2_s > 0.0);
}

#[test]
fn pressure_mode_requires_endpoint_source_support_and_preserves_atomic_refusal() {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "pressure".into(), abscissa_dims: Pressure::DIMS,
            knots: vec![(95000.0, 2e-5), (105000.0, 2e-5)], dims: BINARY_DIFFUSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with_quantity("pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            95000.0, 105000.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "bounded pressure fixture".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    let narrow = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic pair".into(),
        phase: "gas".into(), process: "pressure fixture".into(), revision: 0 }, claims, vec![]).unwrap();
    let context = point(400.0, 100000.0);
    let (a, b) = (gas_at(0.5, 120000.0), gas_at(0.5, 80000.0));
    let saved = (a.clone(), b.clone()); let pairs = [pair("A", "B", &narrow, &context)];
    // The mean face is supported, but neither physical endpoint is.
    a.diffusion_flux_to(&b, geometry(), &pairs, pressure_controls(), || false).unwrap();
    assert!(a.diffuse_with_driving(&b, geometry(), 0.01, &pairs,
        pressure_controls(), driving(), || false).is_err());
    let d = card(2e-5); let pairs = [pair("A", "B", &d, &context)];
    let mut calls = 0;
    a.diffuse_with_driving(&b, geometry(), 0.01, &pairs,
        pressure_controls(), driving(), || { calls += 1; false }).unwrap();
    for stop in [1, calls/2, calls] {
        let mut count = 0;
        assert_eq!(a.diffuse_with_driving(&b, geometry(), 0.01, &pairs,
            pressure_controls(), driving(), || { count += 1; count == stop }),
            Err(GasEnergyError::Cancelled));
    }
    assert!(a.diffuse_with_driving(&b, geometry(), 1e6, &pairs,
        pressure_controls(), driving(), || false).is_err());
    assert!(a.diffusion_flux_to_with_driving(&b, geometry(), &pairs, pressure_controls(),
        GasDiffusionDriving::CompositionAndPressureMass { mass_flux_tolerance_kg_m2_s: f64::NAN },
        || false).is_err());
    assert_eq!((a, b), saved);
}

#[test]
fn pressure_diffusion_pair_fixed_network_and_adaptive_window_conserve_and_increase_entropy() {
    use molecular::adaptive::{advance_candidate_with_driving, DiffusionNetworkAdvanceControls};
    let gases = BTreeMap::from([(0, gas_at(0.5, 120000.0)), (1, gas_at(0.5, 80000.0))]);
    let saved = gases.clone(); let c = pressure_controls();
    let endpoints = gases.keys().map(|&id| (id, c.exchange.gas[0])).collect();
    let solve = molecular::DiffusionNetworkSolve { exchange: SpeciesNetworkSolve {
        maximum_species_fraction: c.exchange.maximum_species_fraction,
        mass_tolerance_kg: c.exchange.mass_tolerance_kg, energy_tolerance_j: c.exchange.energy_tolerance_j },
        entropy_tolerance_j_k: c.entropy_tolerance_j_k };
    let d = card(2e-5); let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let links = [molecular::DiffusionNetworkLink { id: 1, parcels: [0, 1], geometry: geometry(),
        pairs: &pairs, solve: molecular::DiffusionFaceSolve {
            maximum_temperature_difference_k: c.maximum_temperature_difference_k,
            maximum_pressure_difference_pa: c.maximum_pressure_difference_pa,
            gradient_tolerance_per_m: c.gradient_tolerance_per_m,
            molar_flux_tolerance_mol_m2_s: c.molar_flux_tolerance_mol_m2_s } }];
    let pair = gases[&0].diffuse_with_driving(&gases[&1], geometry(), 0.01, &pairs,
        c, driving(), || false).unwrap();
    let graph = molecular::step_candidate_with_driving(&gases, &links, 0.01,
        &endpoints, solve, driving(), || false).unwrap();
    assert!(pair.entropy_change_j_k() > 0.0);
    assert!(graph.step().entropy_change_j_k > 0.0);
    assert!(graph.step().exchange.maximum_species_mass_residual_kg <= c.exchange.mass_tolerance_kg);
    assert!(graph.step().exchange.total_energy_residual_j.abs() <= c.exchange.energy_tolerance_j);
    for (id, expected) in [(0, pair.exchange().left()), (1, pair.exchange().right())] {
        let actual = &graph.parcels()[&id];
        near(actual.mass_kg(), gases[&id].mass_kg(), c.exchange.mass_tolerance_kg);
        near(actual.state().state().temperature, expected.state().state().temperature, 1e-7);
        near(actual.state().state().pressure, expected.state().state().pressure, 1e-5);
    }
    let mut reverse = links.clone(); reverse[0].parcels = [1, 0];
    let reversed = molecular::step_candidate_with_driving(&gases, &reverse, 0.01,
        &endpoints, solve, driving(), || false).unwrap();
    assert_eq!(graph.parcels(), reversed.parcels());
    let window = advance_candidate_with_driving(&gases, &links, 0.1, &endpoints, solve,
        DiffusionNetworkAdvanceControls { maximum_step_s: 0.05, minimum_step_s: 1e-6,
            maximum_evaluations: 1000, absolute_species_tolerance_kg: 1e-9,
            relative_species_tolerance: 1e-5, temperature_tolerance_k: 1e-3 },
        driving(), || false).unwrap();
    let report = window.report();
    assert_eq!(report.duration_s, 0.1);
    assert!(report.accepted_substeps >= 4);
    assert!(report.entropy_change_j_k > 0.0);
    assert!(report.maximum_species_mass_residual_kg <= c.exchange.mass_tolerance_kg);
    assert!(report.total_energy_residual_j.abs() <= c.exchange.energy_tolerance_j);
    assert!(report.endpoint_energy_change_j.abs() <= c.exchange.energy_tolerance_j);
    assert_eq!(report.last_substep.links[0].flux.driving, driving());
    for (&id, parcel) in window.parcels() {
        near(parcel.mass_kg(), gases[&id].mass_kg(), c.exchange.mass_tolerance_kg);
    }
    assert_eq!(gases, saved);
}
