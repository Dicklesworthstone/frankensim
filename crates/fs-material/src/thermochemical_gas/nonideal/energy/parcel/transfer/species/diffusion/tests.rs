//! G1/G3 source-bound real-gas diffusion and finite original-owner balances.
mod adaptive;
mod thermal;

use super::super::tests::{at, close, energy, real_model, source, species_mass};
use super::*;
use crate::state_point::MaterialPropertySelection;
use crate::thermochemical_gas::tests::point;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use network::{RealGasDiffusionNetworkLink, RealGasDiffusionNetworkSolve};
use std::collections::BTreeMap;

fn controls() -> RealGasDiffusionSolve {
    RealGasDiffusionSolve {
        maximum_temperature_difference_k: 1.0,
        maximum_pressure_difference_pa: 1e6,
        diffusion: MaxwellStefanSolveV1::default(),
        entropy_pairing_tolerance_w_m3_k: 1e-10,
    }
}
fn physical() -> RealGasDiffusionNetworkSolve {
    RealGasDiffusionNetworkSolve {
        exchange: super::super::network::SpeciesNetworkSolve {
            maximum_species_fraction: 0.2,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-7,
        },
        entropy_tolerance_j_k: 1e-9,
    }
}
fn geometry() -> GasDiffusionGeometry {
    GasDiffusionGeometry {
        area_m2: 1e-4,
        distance_m: 0.02,
    }
}
fn pair(references: [f64; 2]) -> [IsochoricRealGasParcel; 2] {
    [
        at(real_model(source(0.75, references), 0.03), 330.0, 4e6, 11.0),
        at(real_model(source(0.35, references), 0.03), 330.0, 4e6, 22.0),
    ]
}
fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 3e-8), (450.0, 3e-8)],
                dims: BINARY_DIFFUSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    250.0,
                    450.0,
                )
                .with_quantity(
                    "pressure",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::Pressure,
                        ValueForm::Static,
                    )),
                    1000.0,
                    5e6,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored dense-gas MS test coefficient".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic A/B pair".into(),
            phase: "gas".into(),
            process: "real-gas diffusion fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn binding<'a>(card: &'a MaterialCard, context: &'a QueryPoint) -> BinaryDiffusionMaterial<'a> {
    BinaryDiffusionMaterial {
        species: [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        material: card,
        context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn endpoints(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
) -> BTreeMap<u64, super::super::GasEnergySolve> {
    parcels.keys().map(|&id| (id, energy())).collect()
}

#[test]
fn source_bound_face_matches_independent_pr_binary_values_and_reverses() {
    let [left, right] = pair([0.0; 2]);
    let card = card();
    let context = point(330.0, 4e6);
    let model =
        RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false).unwrap();
    let mut work = model.workspace().unwrap();
    let flux = model
        .face_flux(&left, &right, geometry(), controls(), &mut work, || false)
        .unwrap();
    close(flux.molar_concentration_mol_m3, 1514.9804470555855, 1e-9);
    close(
        flux.species[0].molar_flux_mol_m2_s,
        0.0002253251762677151,
        1e-14,
    );
    close(
        flux.force_entropy_production_w_m3_k,
        0.6938842813864166,
        1e-10,
    );
    assert!(flux.drag_entropy_production_w_m3_k > 0.0);
    close(
        flux.force_entropy_production_w_m3_k,
        flux.drag_entropy_production_w_m3_k,
        1e-10,
    );
    assert_eq!(flux.pairs[0].sample.value_si, 3e-8);
    assert!(
        flux.pairs[0]
            .sample
            .face
            .query_point()
            .contains(&("pressure".into(), 4e6))
    );
    let reverse = model
        .face_flux(&right, &left, geometry(), controls(), &mut work, || false)
        .unwrap();
    for i in 0..2 {
        close(
            reverse.species[i].molar_flux_mol_m2_s,
            -flux.species[i].molar_flux_mol_m2_s,
            1e-14,
        );
    }
    close(
        reverse.force_entropy_production_w_m3_k,
        flux.force_entropy_production_w_m3_k,
        1e-10,
    );
}

#[test]
fn finite_network_preserves_original_species_energy_and_actual_entropy() {
    let [left, right] = pair([0.0; 2]);
    let card = card();
    let context = point(330.0, 4e6);
    let model =
        RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false).unwrap();
    let links = [RealGasDiffusionNetworkLink {
        id: 7,
        parcels: [10, 20],
        model: &model,
        geometry: geometry(),
        solve: controls(),
    }];
    let parcels = BTreeMap::from([(10, left), (20, right)]);
    let saved = parcels.clone();
    let candidate = network::step_candidate(
        &parcels,
        &links,
        1.0,
        &endpoints(&parcels),
        physical(),
        || false,
    )
    .unwrap();
    let flux = candidate.step().links[0].flux.as_ref().unwrap();
    for i in 0..2 {
        let row = &flux.species[i];
        close(
            species_mass(&candidate.parcels()[&10], &row.species)
                - species_mass(&parcels[&10], &row.species),
            -row.mass_rate_left_to_right_kg_s,
            1e-14,
        );
        close(
            species_mass(&candidate.parcels()[&20], &row.species)
                - species_mass(&parcels[&20], &row.species),
            row.mass_rate_left_to_right_kg_s,
            1e-14,
        );
    }
    let mut energy_change = 0.0;
    let mut entropy_change = 0.0;
    for (&id, old) in &parcels {
        let new = &candidate.parcels()[&id];
        energy_change += new.mass_kg() * new.state().specific_internal_energy_j_kg()
            - old.mass_kg() * old.state().specific_internal_energy_j_kg();
        entropy_change += new.mass_kg() * new.state().specific_entropy_j_kg_k().unwrap()
            - old.mass_kg() * old.state().specific_entropy_j_kg_k().unwrap();
        assert_eq!(old.volume_m3(), new.volume_m3());
        close(new.mass_kg(), old.mass_kg(), 1e-14);
    }
    close(energy_change, 0.0, 1e-7);
    close(candidate.step().energy_change_j, energy_change, 1e-7);
    close(candidate.step().entropy_change_j_k, entropy_change, 1e-9);
    assert!(candidate.step().entropy_change_j_k > 0.0);
    assert_eq!(parcels, saved);
}

#[test]
fn parallel_face_order_and_formation_datums_do_not_change_physical_transport() {
    let card = card();
    let context = point(330.0, 4e6);
    let mut outputs = Vec::new();
    for references in [[0.0; 2], [1e19, -2e19]] {
        let [left, right] = pair(references);
        let model = RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false)
            .unwrap();
        let mut links = [
            RealGasDiffusionNetworkLink {
                id: 7,
                parcels: [10, 20],
                model: &model,
                geometry: geometry(),
                solve: controls(),
            },
            RealGasDiffusionNetworkLink {
                id: 4,
                parcels: [10, 20],
                model: &model,
                geometry: geometry(),
                solve: controls(),
            },
        ];
        let parcels = BTreeMap::from([(10, left), (20, right)]);
        let first = network::step_candidate(
            &parcels,
            &links,
            0.5,
            &endpoints(&parcels),
            physical(),
            || false,
        )
        .unwrap();
        links.reverse();
        let second = network::step_candidate(
            &parcels,
            &links,
            0.5,
            &endpoints(&parcels),
            physical(),
            || false,
        )
        .unwrap();
        assert_eq!(first, second);
        outputs.push(first);
    }
    for id in [10, 20] {
        close(
            outputs[0].parcels()[&id].state().temperature_k(),
            outputs[1].parcels()[&id].state().temperature_k(),
            1e-10,
        );
        close(
            outputs[0].parcels()[&id].state().pressure_pa(),
            outputs[1].parcels()[&id].state().pressure_pa(),
            1e-5,
        );
    }
    close(
        outputs[0].step().entropy_change_j_k,
        outputs[1].step().entropy_change_j_k,
        1e-10,
    );
}

#[test]
fn equilibrium_closed_faces_and_source_refusals_are_explicit() {
    let [left, right] = pair([0.0; 2]);
    let card = card();
    let context = point(330.0, 4e6);
    let model =
        RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false).unwrap();
    let mut work = model.workspace().unwrap();
    let zero = model
        .face_flux(&left, &left, geometry(), controls(), &mut work, || false)
        .unwrap();
    assert!(
        zero.species
            .iter()
            .all(|row| row.molar_flux_mol_m2_s == 0.0)
    );
    assert_eq!(zero.force_entropy_production_w_m3_k, 0.0);
    let changed = at(real_model(source(0.35, [0.0; 2]), 0.04), 330.0, 4e6, 22.0);
    assert!(
        model
            .face_flux(&left, &changed, geometry(), controls(), &mut work, || false)
            .is_err()
    );
    let hot = at(real_model(source(0.35, [0.0; 2]), 0.03), 335.0, 4e6, 22.0);
    assert!(
        model
            .face_flux(&left, &hot, geometry(), controls(), &mut work, || false)
            .is_err()
    );
    let pure = left
        .model()
        .source()
        .with_mole_amounts(&[fs_qty::Amount::new(1.0), fs_qty::Amount::new(0.0)])
        .unwrap();
    let one_sided = at(real_model(pure, 0.03), 330.0, 4e6, 22.0);
    assert!(
        model
            .face_flux(&left, &one_sided, geometry(), controls(), &mut work, || {
                false
            })
            .is_err()
    );
    let outside = at(real_model(source(0.35, [0.0; 2]), 0.03), 330.0, 8e6, 22.0);
    let parcels = BTreeMap::from([(10, left), (20, outside)]);
    let links = [RealGasDiffusionNetworkLink {
        id: 1,
        parcels: [10, 20],
        model: &model,
        geometry: GasDiffusionGeometry {
            area_m2: 0.0,
            ..geometry()
        },
        solve: controls(),
    }];
    let candidate = network::step_candidate(
        &parcels,
        &links,
        1.0,
        &endpoints(&parcels),
        physical(),
        || false,
    )
    .unwrap();
    assert_eq!(candidate.parcels(), &parcels);
    assert!(candidate.step().links[0].flux.is_none());
    let _ = right;
}

#[test]
fn late_cancellation_and_gross_withdrawal_return_no_candidate() {
    let [left, right] = pair([0.0; 2]);
    let card = card();
    let context = point(330.0, 4e6);
    let model =
        RealGasDiffusionModel::new(left.model(), &[binding(&card, &context)], || false).unwrap();
    let links = [RealGasDiffusionNetworkLink {
        id: 1,
        parcels: [10, 20],
        model: &model,
        geometry: geometry(),
        solve: controls(),
    }];
    let parcels = BTreeMap::from([(10, left), (20, right)]);
    let saved = parcels.clone();
    let mut polls = 0usize;
    network::step_candidate(
        &parcels,
        &links,
        1.0,
        &endpoints(&parcels),
        physical(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0usize;
    assert!(matches!(
        network::step_candidate(
            &parcels,
            &links,
            1.0,
            &endpoints(&parcels),
            physical(),
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(network::RealGasDiffusionNetworkError::Cancelled)
    ));
    assert!(
        network::step_candidate(
            &parcels,
            &links,
            1e9,
            &endpoints(&parcels),
            physical(),
            || false
        )
        .is_err()
    );
    assert_eq!(parcels, saved);
}

#[test]
fn the_face_audit_refuses_a_net_twice_its_bound_and_admits_the_published_residual() {
    let u = f64::EPSILON / 2.0;
    let (area, duration) = (1e-4, 20.0);
    // With no published residual the bound is 5u sum|m_i|; the allowance adds
    // 6u sum|m_i| for the network's three compensated sums.
    let gross = 2e-6;
    let bound = 5.0 * u * gross;
    let balanced = network::face_allowance(&[1e-6, -1e-6], 0.0, area, duration).unwrap();
    assert_eq!(balanced, bound + 6.0 * u * gross);
    // Planted negative (a): a face netting about twice its bound is refused.
    let over = [1e-6, -(1e-6 - 2.0 * bound)];
    let ratio = (over[0] + over[1]) / bound;
    assert!((1.8..2.2).contains(&ratio), "{ratio}");
    assert!(network::face_allowance(&over, 0.0, area, duration).is_err());
    // The owner's residual widens the bound by |r| A dt, and no further.
    let residual = 1e-14;
    let rate = residual * area * duration;
    let within = [1e-6, -(1e-6 - rate)];
    let beyond = [1e-6, -(1e-6 - 2.0 * rate)];
    let allowance = network::face_allowance(&within, residual, area, duration).unwrap();
    assert!(allowance >= (within[0] + within[1]).abs());
    assert!(network::face_allowance(&within, -residual, area, duration).is_ok());
    assert!(network::face_allowance(&beyond, residual, area, duration).is_err());
    assert!(network::face_allowance(&within, f64::NAN, area, duration).is_err());
}

#[test]
fn a_general_path_mass_frame_face_is_admitted_within_its_published_residual() {
    use fs_qty::{AmountConcentration, MolarMass};
    use fs_thermochem::transport::maxwell_stefan::MaxwellStefanV1;
    // Planted case (b1). Three species, all present, with three different pair
    // diffusivities take the owner's dense LU path, where the mass frame is one
    // row of the solve rather than an explicit closure (fs-thermochem
    // transport/maxwell_stefan/composition.rs): n = 3 passes the two-species
    // closure (:23), unequal pairs clear the common diffusivity (:50, :58), and
    // three positive fractions pass the one-carrier closure (:107) to the LU
    // solve (:159). Its published residual is that row's rounding.
    let owner = MaxwellStefanV1::new(vec![
        (SpeciesId::new("A").unwrap(), MolarMass::new(0.028)),
        (SpeciesId::new("B").unwrap(), MolarMass::new(0.004)),
        (SpeciesId::new("C").unwrap(), MolarMass::new(0.044)),
    ])
    .unwrap();
    let mut work = owner.workspace().unwrap();
    let d = [0.0, 3e-8, 1e-8, 3e-8, 0.0, 7e-8, 1e-8, 7e-8, 0.0];
    let (area, duration) = (1e-4, 20.0);
    let mut nonzero = 0;
    for x in [
        [0.5, 0.3, 0.2],
        [0.45, 0.35, 0.2],
        [0.6, 0.15, 0.25],
        [0.3, 0.3, 0.4],
    ] {
        for forces in [[-37.0, 48.0, -11.0], [5.0, -3.0, -2.0], [-1.5, -2.5, 4.0]] {
            let flux = owner
                .evaluate_driving_forces(
                    &x,
                    &forces,
                    AmountConcentration::new(1500.0),
                    &d,
                    MaxwellStefanSolveV1::default(),
                    &mut work,
                    || false,
                )
                .unwrap();
            let residual = flux.mass_frame_residual_kg_m2_s();
            let scale: f64 = flux.mass_fluxes().iter().map(|j| j.abs()).sum();
            assert!(residual.abs() <= 1e-16 + 1e-10 * scale);
            if residual == 0.0 {
                continue;
            }
            nonzero += 1;
            let masses: Vec<_> = flux
                .mass_fluxes()
                .iter()
                .map(|j| j * area * duration)
                .collect();
            let net = masses
                .iter()
                .fold(fs_math::dd::Dd::ZERO, |sum, m| {
                    sum + fs_math::dd::Dd::from_f64(*m)
                })
                .to_f64();
            let allowance = network::face_allowance(&masses, residual, area, duration).unwrap();
            assert!(allowance >= net.abs(), "{allowance} < {net}");
        }
    }
    assert!(
        nonzero > 0,
        "no general-path case published a nonzero residual"
    );
}
