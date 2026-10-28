//! Synthetic analytical limits through actual material and gas owners.
use super::*;
use crate::porous::permeability::{
    DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
};
use crate::state_point::MaterialPropertySelection;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    tests::{model, point, species, transport},
};
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, network::snapshot::SnapshotNozzleLink,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::diffusion::{
    BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkError, SpeciesNetworkSolve,
    transport::{
        GasTransportError, GasTransportPaths, GasTransportSolve, step_candidate_with_paths,
        dusty_gas::DustyGasNetworkLink,
    },
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{
    Dims, Pressure, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    mixture::FrozenIdealGasMixtureModelV1,
};
use std::collections::BTreeMap;

const T: f64 = 400.0;
const P: f64 = 100_000.0;
const LENGTH: f64 = 0.01;
const AREA: f64 = 0.01;
const RADIUS: f64 = 1e-7;
const POROSITY: f64 = 0.4;
const TORTUOSITY: f64 = 2.0;

fn controls() -> DustyGasFaceSolve {
    DustyGasFaceSolve {
        maximum_temperature_difference_k: 1e-8,
        maximum_relative_pressure_difference: 0.25,
        gradient_tolerance_per_m: 1e-8,
    }
}
fn support() -> PorousGasSupport {
    PorousGasSupport {
        temperature_k: (200.0, 1800.0),
        pressure_pa: (1000.0, 1e7),
        maximum_superficial_speed_m_s: 100.0,
    }
}
fn medium(k: f64, length: f64) -> DustyGasMedium {
    let layer = PorousGasLayer::try_new(
        IntrinsicPermeability::isotropic(k).unwrap(),
        length,
        0.0,
        "synthetic homogeneous capillary geometry",
    ).unwrap();
    let barrier = DarcyGasBarrier::try_new(vec![layer], [1.0, 0.0, 0.0], support()).unwrap();
    DustyGasMedium::try_new(
        barrier,
        DustyGasPores::try_new(RADIUS, POROSITY, TORTUOSITY, "synthetic open pore geometry").unwrap(),
    ).unwrap()
}
fn binary(x: f64, pressure: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], 0.0), point(T, pressure), 0.001).unwrap()
}
fn pure(name: &str, pressure: f64) -> ClosedGasParcel {
    let (mass, cp, mu, conductivity) = match name {
        "A" => (0.028, 3.5, 1.7e-5, 0.025),
        "B" => (0.004, 2.5, 2e-5, 0.15),
        _ => unreachable!(),
    };
    let gas = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species(name, mass, cp, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        ).unwrap(),
        vec![transport(name, mu, conductivity)],
        None,
    ).unwrap();
    ClosedGasParcel::new(gas, point(T, pressure), 0.001).unwrap()
}
fn ternary(fractions: Vec<f64>) -> ClosedGasParcel {
    let gas = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species("A", 0.028, 3.5, 0.0), species("B", 0.004, 2.5, 0.0),
                species("C", 0.044, 4.0, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        ).unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15),
            transport("C", 1.5e-5, 0.02)],
        None,
    ).unwrap();
    ClosedGasParcel::new(gas, point(T, P), 0.001).unwrap()
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn source(value: PropertyValue, dims: Dims, policy: InterpolationPolicy, pressures: [f64; 2])
    -> MaterialCard
{
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, dims),
        value,
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 200.0, 1800.0)
            .with_quantity("pressure", typed(QuantityKind::Pressure), pressures[0], pressures[1]),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: policy,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic binary gas diffusivity".into(), license: "CC0-1.0".into(), artifact: None,
        },
    }).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic species pair".into(), phase: "gas".into(),
            process: "dusty gas analytic fixture".into(), revision: 0,
        },
        claims, vec![],
    ).unwrap()
}
fn card(value: f64) -> MaterialCard {
    source(
        PropertyValue::Scalar { value, dims: BINARY_DIFFUSIVITY_DIMS },
        BINARY_DIFFUSIVITY_DIMS, InterpolationPolicy::ConstantWithinValidity, [1000.0, 1e7],
    )
}
fn pair<'a>(a: &str, b: &str, card: &'a MaterialCard, context: &'a QueryPoint)
    -> BinaryDiffusionMaterial<'a>
{
    BinaryDiffusionMaterial {
        species: [SpeciesId::new(a).unwrap(), SpeciesId::new(b).unwrap()],
        material: card, context, selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn near(actual: f64, expected: f64) {
    assert!((actual-expected).abs() <= 3e-11*expected.abs().max(1e-12),
        "{actual:.17e} != {expected:.17e}");
}
fn knudsen(molar_mass: f64) -> f64 {
    POROSITY/TORTUOSITY * (2.0*RADIUS/3.0)
        * (8.0*UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*T/(core::f64::consts::PI*molar_mass)).sqrt()
}

#[test]
fn pure_species_recovers_knudsen_plus_pressure_squared_darcy_and_thickness_scaling() {
    let (left, right) = (pure("A", 110_000.0), pure("A", 90_000.0));
    let context = point(300.0, 80_000.0);
    let k = 2e-14;
    let slab = medium(k, LENGTH);
    let flux = left.dusty_gas_flux_to(&right, &slab, AREA, &[], &context, controls(), || false).unwrap();
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
    let gradient = (90_000.0-110_000.0)/LENGTH;
    let expected_diffusion = -knudsen(0.028)*gradient/(r*T);
    // The existing pure transport fixture doubles between 200 K and 2000 K.
    let viscosity = 1.7e-5 * (1.0 + (T - 200.0) / 1800.0);
    let expected_viscous = -k*P*gradient/(viscosity*r*T);
    assert_eq!(flux.species.len(), 1);
    near(flux.species[0].knudsen_diffusivity_m2_s, knudsen(0.028));
    near(flux.species[0].diffusive_molar_flux_mol_m2_s, expected_diffusion);
    near(flux.species[0].viscous_molar_flux_mol_m2_s, expected_viscous);
    near(flux.species[0].molar_flux_mol_m2_s, expected_diffusion+expected_viscous);
    near(flux.species[0].mass_rate_left_to_right_kg_s,
        (expected_diffusion+expected_viscous)*0.028*AREA);
    near(flux.net_molar_flux_mol_m2_s, expected_diffusion+expected_viscous);
    assert!(flux.local_entropy_production_w_m3_k > 0.0);
    assert!(flux.maximum_gradient_residual_per_m < controls().gradient_tolerance_per_m);
    let reverse = right.dusty_gas_flux_to(&left, &slab, AREA, &[], &context, controls(), || false).unwrap();
    near(reverse.species[0].molar_flux_mol_m2_s, -flux.species[0].molar_flux_mol_m2_s);
    near(reverse.local_entropy_production_w_m3_k, flux.local_entropy_production_w_m3_k);
    let thick = left.dusty_gas_flux_to(&right, &medium(k, 2.0*LENGTH), AREA,
        &[], &context, controls(), || false).unwrap();
    near(thick.species[0].mass_rate_left_to_right_kg_s,
        0.5*flux.species[0].mass_rate_left_to_right_kg_s);
    let zero_darcy = left.dusty_gas_flux_to(&right, &medium(0.0, LENGTH), AREA,
        &[], &context, controls(), || false).unwrap();
    near(zero_darcy.species[0].molar_flux_mol_m2_s, expected_diffusion);
    assert_eq!(zero_darcy.species[0].viscous_molar_flux_mol_m2_s, 0.0);
}

#[test]
fn binary_isobaric_graham_counterdiffusion_has_nonzero_total_molar_flux() {
    let (left, right) = (binary(0.75, P), binary(0.25, P));
    let (d, context, slab) = (card(2e-5), point(T, P), medium(0.0, LENGTH));
    let bindings = [pair("A", "B", &d, &context)];
    let flux = left.dusty_gas_flux_to(&right, &slab, AREA, &bindings, &context, controls(), || false).unwrap();
    let concentration = P/(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*T);
    let diffusion = 2e-5*POROSITY/TORTUOSITY;
    let (ka, kb) = (knudsen(0.028), knudsen(0.004));
    let denominator = diffusion + 0.5*ka + 0.5*kb;
    let gradient = (0.25-0.75)/LENGTH;
    let na = -concentration*ka*diffusion*gradient/denominator;
    let nb = concentration*kb*diffusion*gradient/denominator;
    near(flux.species[0].molar_flux_mol_m2_s, na);
    near(flux.species[1].molar_flux_mol_m2_s, nb);
    near(na/nb, -(0.004_f64/0.028).sqrt());
    near(flux.net_molar_flux_mol_m2_s, na+nb);
    assert!(flux.net_molar_flux_mol_m2_s.abs() > 1e-4);
    near(flux.pairs[0].diffusivity_m2_s, 2e-5);
    for s in &flux.species {
        assert_eq!(s.viscous_molar_flux_mol_m2_s, 0.0);
        near(s.face_mole_fraction, 0.5);
    }
    let reversed = [pair("B", "A", &d, &context)];
    let backward = right.dusty_gas_flux_to(&left, &slab, AREA, &reversed,
        &context, controls(), || false).unwrap();
    for i in 0..2 { near(backward.species[i].molar_flux_mol_m2_s, -flux.species[i].molar_flux_mol_m2_s); }
}

#[test]
fn ternary_cross_friction_can_drive_a_species_up_its_own_partial_pressure_gradient() {
    let (left, right) = (ternary(vec![0.6, 0.2, 0.2]), ternary(vec![0.2, 0.6005, 0.1995]));
    let cards = [card(1e-5), card(2e-5), card(5e-5)];
    let context = point(T, P);
    let slab = medium(0.0, LENGTH);
    let bindings = [pair("A", "B", &cards[0], &context), pair("A", "C", &cards[1], &context),
        pair("B", "C", &cards[2], &context)];
    let flux = left.dusty_gas_flux_to(&right, &slab, AREA, &bindings, &context, controls(), || false).unwrap();
    // Independent Cramer determinants evaluated with 60-digit decimal arithmetic;
    // no copied production LU, gauge row, or binary effective-D substitution.
    let expected = [0.0013174171040320678, -0.003465053821027934, -0.000006182251680518990];
    for (actual, expected) in flux.species.iter().zip(expected) {
        near(actual.molar_flux_mol_m2_s, expected);
    }
    assert!(flux.species[2].molar_flux_mol_m2_s < 0.0);
    // C moves right -> left, into its higher partial pressure.
    assert!(right.model().thermodynamics().mole_composition().fractions()[2]
        < left.model().thermodynamics().mole_composition().fractions()[2]);
    assert!(flux.local_entropy_production_w_m3_k > 0.0);
    let reordered = [bindings[2].clone(), bindings[0].clone(), bindings[1].clone()];
    let again = left.dusty_gas_flux_to(&right, &slab, AREA, &reordered, &context, controls(), || false).unwrap();
    for (a, b) in flux.species.iter().zip(&again.species) {
        assert_eq!(a.species, b.species);
        assert_eq!(a.molar_flux_mol_m2_s, b.molar_flux_mol_m2_s);
    }
}

#[test]
fn face_composition_uses_concentrations_and_pressure_dependent_binary_sources_are_allowed() {
    let (left, right) = (binary(0.8, 110_000.0), binary(0.2, 90_000.0));
    let context = point(T, P);
    let d = source(PropertyValue::Curve {
        abscissa: "pressure".into(), abscissa_dims: Pressure::DIMS,
        knots: vec![(80_000.0, 2.5e-5), (120_000.0, 1.5e-5)],
        dims: BINARY_DIFFUSIVITY_DIMS,
    }, BINARY_DIFFUSIVITY_DIMS, InterpolationPolicy::LinearInside, [80_000.0, 120_000.0]);
    let bindings = [pair("A", "B", &d, &context)];
    let slab = medium(2e-14, LENGTH);
    let flux = left.dusty_gas_flux_to(&right, &slab, AREA, &bindings, &context, controls(), || false).unwrap();
    near(flux.species[0].face_mole_fraction, (110_000.0*0.8+90_000.0*0.2)/200_000.0);
    near(flux.species[1].face_mole_fraction, 0.47);
    near(flux.pairs[0].diffusivity_m2_s, 2e-5);
    assert!(flux.face_state.is_some());
    assert!(flux.darcy.is_some());
}

#[test]
fn source_endpoint_support_pair_identity_closed_faces_and_cancellation_are_explicit() {
    let (left, right) = (binary(0.75, 110_000.0), binary(0.25, 90_000.0));
    let saved = (left.clone(), right.clone());
    let (d, context, slab) = (card(2e-5), point(T, P), medium(0.0, LENGTH));
    let narrow = source(PropertyValue::Scalar { value: 2e-5, dims: BINARY_DIFFUSIVITY_DIMS },
        BINARY_DIFFUSIVITY_DIMS, InterpolationPolicy::ConstantWithinValidity, [95_000.0, 105_000.0]);
    assert!(left.dusty_gas_flux_to(&right, &slab, AREA,
        &[pair("A", "B", &narrow, &context)], &context, controls(), || false).is_err());
    let wrong_dims = Dims([2, 0, 0, 0, 0, 0]);
    let wrong = source(PropertyValue::Scalar { value: 2e-5, dims: wrong_dims },
        wrong_dims, InterpolationPolicy::ConstantWithinValidity, [1000.0, 1e7]);
    assert!(left.dusty_gas_flux_to(&right, &slab, AREA,
        &[pair("A", "B", &wrong, &context)], &context, controls(), || false).is_err());
    assert!(left.dusty_gas_flux_to(&right, &slab, AREA, &[], &context, controls(), || false).is_err());
    assert!(left.dusty_gas_flux_to(&right, &slab, AREA,
        &[pair("A", "C", &d, &context)], &context, controls(), || false).is_err());
    let closed = left.dusty_gas_flux_to(&right, &slab, 0.0, &[], &context, controls(), || false).unwrap();
    assert!(closed.species.is_empty());
    assert!(closed.pairs.is_empty());
    assert!(closed.face_state.is_none());
    assert!(closed.darcy.is_none());
    assert!(matches!(left.dusty_gas_flux_to(&right, &slab, AREA,
        &[pair("A", "B", &d, &context)], &context, controls(), || true),
        Err(GasEnergyError::Cancelled)));
    let strict = DustyGasFaceSolve { maximum_relative_pressure_difference: 0.01, ..controls() };
    assert!(left.dusty_gas_flux_to(&right, &slab, AREA,
        &[pair("A", "B", &d, &context)], &context, strict, || false).is_err());
    assert_eq!((&left, &right), (&saved.0, &saved.1));
}

#[test]
fn pure_opposite_species_need_no_artificial_trace_fraction_and_equilibrium_has_zero_flux() {
    let (left, right) = (pure("A", P), pure("B", P));
    let (d, context, slab) = (card(2e-5), point(T, P), medium(0.0, LENGTH));
    let bindings = [pair("A", "B", &d, &context)];
    let flux = left.dusty_gas_flux_to(&right, &slab, AREA, &bindings, &context, controls(), || false).unwrap();
    assert_eq!(flux.species.len(), 2);
    assert!(flux.species[0].molar_flux_mol_m2_s > 0.0);
    assert!(flux.species[1].molar_flux_mol_m2_s < 0.0);
    near(flux.species[0].molar_flux_mol_m2_s/flux.species[1].molar_flux_mol_m2_s,
        -(0.004_f64/0.028).sqrt());
    let same = binary(0.75, P);
    let equilibrium = same.dusty_gas_flux_to(&same, &slab, AREA, &bindings, &context, controls(), || false).unwrap();
    for s in equilibrium.species { assert_eq!(s.molar_flux_mol_m2_s, 0.0); }
    assert_eq!(equilibrium.local_entropy_production_w_m3_k, 0.0);
}

fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn network_controls() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1, mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    }
}
fn link<'a>(medium: &'a DustyGasMedium, pairs: &'a [BinaryDiffusionMaterial<'a>], context: &'a QueryPoint)
    -> DustyGasNetworkLink<'a>
{
    DustyGasNetworkLink {
        id: 91, parcels: [0, 1], medium, area_m2: AREA, pairs,
        gas_context: context, solve: controls(),
    }
}

#[test]
fn actual_combined_owner_conserves_species_and_donor_partial_enthalpy_with_unequal_molar_flow() {
    let parcels = BTreeMap::from([(0, pure("A", P)), (1, pure("B", P))]);
    let saved = parcels.clone();
    let endpoints = parcels.keys().map(|&id| (id, energy())).collect();
    let (d, context, slab) = (card(2e-5), point(T, P), medium(0.0, LENGTH));
    let bindings = [pair("A", "B", &d, &context)];
    let links = [link(&slab, &bindings, &context)];
    let candidate = step_candidate_with_paths(
        &parcels, GasTransportPaths { dusty_gas: &links, ..Default::default() },
        0.01, &endpoints, network_controls(), || false,
    ).unwrap();
    let step = candidate.step();
    assert_eq!(step.exchange.fluxes.len(), 2);
    assert!(step.exchange.total_mass_residual_kg.abs() < 1e-13);
    assert!(step.exchange.maximum_species_mass_residual_kg < 1e-13);
    assert!(step.exchange.total_energy_residual_j.abs() < 1e-8);
    assert!(step.entropy_change_j_k > 0.0);
    for row in &step.exchange.fluxes {
        let donor = &parcels[&row.transfer.donor];
        let species = &donor.model().thermodynamics().components()[0];
        let enthalpy = species.evaluate(fs_qty::Temperature::new(T)).unwrap()
            .properties().h().value()/species.molar_mass().value();
        near(row.source_enthalpy_j_kg, enthalpy);
        near(row.carried_enthalpy_j, row.transfer.mass_kg*enthalpy);
    }
    for next in candidate.parcels().values() {
        assert_eq!(next.model().thermodynamics().components().len(), 2);
    }
    assert_eq!(parcels, saved);
}

#[test]
fn dusty_gas_and_nozzle_share_the_original_species_withdrawal_budget() {
    let parcels = BTreeMap::from([
        (0, binary(0.1, P)), (1, binary(0.01, P)), (2, binary(0.1, 50_000.0)),
    ]);
    let saved = parcels.clone();
    let endpoints = parcels.keys().map(|&id| (id, energy())).collect();
    let (d, context, slab) = (card(2e-5), point(T, P), medium(0.0, LENGTH));
    let bindings = [pair("A", "B", &d, &context)];
    let mut pore = link(&slab, &bindings, &context);
    let dt = 0.01;
    let flux = parcels[&0].dusty_gas_flux_to(&parcels[&1], &slab, AREA, &bindings,
        &context, controls(), || false).unwrap();
    let initial_a = parcels[&0].mass_kg()*0.1*0.028/(0.1*0.028+0.9*0.004);
    pore.area_m2 *= 0.06*initial_a/(dt*flux.species[0].mass_rate_left_to_right_kg_s);
    let mut nozzle = SnapshotNozzleLink {
        id: 7, parcels: [0, 2],
        nozzle: GasNozzle { area_m2: 1e-8, discharge_coefficient: 0.8 },
        nozzle_solve: energy(),
    };
    let rate = parcels[&0].nozzle_flow_to(50_000.0, nozzle.nozzle, energy(), || false).unwrap().mass_flow_kg_s;
    nozzle.nozzle.area_m2 *= 0.06*parcels[&0].mass_kg()/(rate*dt);
    step_candidate_with_paths(&parcels,
        GasTransportPaths { dusty_gas: &[pore.clone()], ..Default::default() },
        dt, &endpoints, network_controls(), || false).unwrap();
    step_candidate_with_paths(&parcels,
        GasTransportPaths { nozzles: &[nozzle], ..Default::default() },
        dt, &endpoints, network_controls(), || false).unwrap();
    assert!(matches!(step_candidate_with_paths(&parcels,
        GasTransportPaths { dusty_gas: &[pore], nozzles: &[nozzle], ..Default::default() },
        dt, &endpoints, network_controls(), || false),
        Err(GasTransportError::Exchange(SpeciesNetworkError::Withdrawal { id: 0, .. }))));
    assert_eq!(parcels, saved);
}

#[test]
fn medium_and_controls_refuse_layered_pores_slip_double_counting_and_nonphysical_inputs() {
    for (radius, porosity, tortuosity) in [
        (0.0, POROSITY, TORTUOSITY),
        (-RADIUS, POROSITY, TORTUOSITY),
        (f64::NAN, POROSITY, TORTUOSITY),
        (f64::INFINITY, POROSITY, TORTUOSITY),
        (RADIUS, 0.0, TORTUOSITY),
        (RADIUS, -0.1, TORTUOSITY),
        (RADIUS, 1.01, TORTUOSITY),
        (RADIUS, f64::NAN, TORTUOSITY),
        (RADIUS, POROSITY, 0.0),
        (RADIUS, POROSITY, 0.5),
        (RADIUS, POROSITY, f64::NAN),
        (RADIUS, POROSITY, f64::INFINITY),
    ] {
        assert!(DustyGasPores::try_new(radius, porosity, tortuosity, "synthetic geometry").is_err());
    }
    assert!(DustyGasPores::try_new(RADIUS, POROSITY, TORTUOSITY, "").is_err());
    let pores = DustyGasPores::try_new(RADIUS, POROSITY, TORTUOSITY, "synthetic geometry").unwrap();
    let layer = |slip| PorousGasLayer::try_new(
        IntrinsicPermeability::isotropic(2e-14).unwrap(), LENGTH, slip,
        "synthetic layer for constitutive admission",
    ).unwrap();
    let layers = DarcyGasBarrier::try_new(vec![layer(0.0), layer(0.0)],
        [1.0, 0.0, 0.0], support()).unwrap();
    // A Darcy series laminate is valid, but its compositions require a resolved
    // multicomponent interface solve before it can claim homogeneous dusty gas.
    assert!(DustyGasMedium::try_new(layers, pores.clone()).is_err());
    let gas = pure("A", P);
    let slip = DarcyGasBarrier::try_new_for_gas_model(vec![layer(5000.0)],
        [1.0, 0.0, 0.0], support(), gas.model()).unwrap();
    // The slip law is independently qualified for this gas; refusing it here
    // prevents adding a second rarefaction correction to Knudsen wall transport.
    assert!(DustyGasMedium::try_new(slip, pores).is_err());

    let slab = medium(2e-14, LENGTH);
    let context = point(T, P);
    let exact = DustyGasFaceSolve {
        maximum_temperature_difference_k: 0.0,
        maximum_relative_pressure_difference: 0.0,
        gradient_tolerance_per_m: 1e-8,
    };
    gas.dusty_gas_flux_to(&gas, &slab, AREA, &[], &context, exact, || false).unwrap();
    for invalid in [
        DustyGasFaceSolve { maximum_temperature_difference_k: -1.0, ..exact },
        DustyGasFaceSolve { maximum_temperature_difference_k: f64::NAN, ..exact },
        DustyGasFaceSolve { maximum_relative_pressure_difference: -0.01, ..exact },
        DustyGasFaceSolve { maximum_relative_pressure_difference: 2.0, ..exact },
        DustyGasFaceSolve { maximum_relative_pressure_difference: f64::INFINITY, ..exact },
        DustyGasFaceSolve { gradient_tolerance_per_m: 0.0, ..exact },
        DustyGasFaceSolve { gradient_tolerance_per_m: f64::NAN, ..exact },
    ] {
        assert!(gas.dusty_gas_flux_to(&gas, &slab, AREA, &[], &context, invalid, || false).is_err());
    }
    for area in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(gas.dusty_gas_flux_to(&gas, &slab, area, &[], &context, exact, || false).is_err());
    }
}

#[test]
fn unrepresentable_nonzero_driving_force_is_refused_instead_of_fabricating_equilibrium() {
    // A zero intrinsic permeability leaves the Knudsen path active and avoids
    // overflowing a Darcy L/k resistance before the diffusion range admission.
    let slab = medium(0.0, 1e308);
    let context = point(T, P);
    let left = pure("A", P);
    let right = pure("A", f64::from_bits(P.to_bits() + 1));
    assert!(right.state().state().pressure > left.state().state().pressure);
    assert!(left.dusty_gas_flux_to(
        &right, &slab, AREA, &[], &context, controls(), || false,
    ).is_err());

    // Equal total pressures can still carry a nonzero composition drive. The
    // exact source mole fractions survive; only their per-length drive cannot
    // be represented. Silently rounding that drive to zero would invent a
    // thermodynamic equilibrium between distinguishable endpoint inventories.
    let left = binary(0.5, P);
    let right = binary(f64::from_bits(0.5_f64.to_bits() + 1), P);
    let xl = left.model().thermodynamics().mole_composition().fractions()[0];
    let xr = right.model().thermodynamics().mole_composition().fractions()[0];
    assert!(xr > xl);
    assert!(right.state().state().pressure*xr > left.state().state().pressure*xl);
    let d = card(2e-5);
    let bindings = [pair("A", "B", &d, &context)];
    assert!(left.dusty_gas_flux_to(
        &right, &slab, AREA, &bindings, &context, controls(), || false,
    ).is_err());
}
