//! Native fixed and adaptive laminate paths; independent synthetic physical controls.
use super::{
    EffusionAdvanceControls, EffusionAdvanceError, advance_candidate,
    advance_candidate_with_row_limit,
};
use crate::porous::permeability::{
    DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
};
use crate::state_point::MaterialPropertySelection;
use crate::thermochemical_gas::energy::parcel::transfer::species::{
    diffusion::{
        BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
        GasDiffusionGeometry,
        dusty_gas::{DustyGasFaceSolve, DustyGasMedium, DustyGasPores, LayeredDustyGasMedium},
        network::{DiffusionFaceSolve, DiffusionNetworkLink},
    },
    inventory,
    network::{
        self, SpeciesNetworkError, SpeciesNetworkSolve, SpeciesNetworkTransfer,
        transport::{
            self as native, GasTransportCandidate, GasTransportError, GasTransportPaths,
            GasTransportSolve,
            dusty_gas::{DustyGasNetworkError, DustyGasNetworkLink, LayeredDustyGasNetworkLink},
        },
    },
};
use crate::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    tests::{model, point, species, species_with_energy_reference, transport},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};
use std::collections::BTreeMap;

// Fixtures only populate the native path field; all sampled physics and state
// reconstruction use the production transport owner below.
fn step_candidate<'a>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &'a [LayeredDustyGasNetworkLink<'a>],
    duration: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasTransportCandidate, GasTransportError> {
    native::step_candidate_with_paths(
        parcels,
        GasTransportPaths {
            layered_dusty_gas: links,
            ..Default::default()
        },
        duration,
        endpoints,
        controls,
        cancelled,
    )
}
#[allow(clippy::too_many_arguments)]
fn step_candidate_with_paths<'a>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    mut paths: GasTransportPaths<'a>,
    links: &'a [LayeredDustyGasNetworkLink<'a>],
    duration: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasTransportCandidate, GasTransportError> {
    assert!(paths.layered_dusty_gas.is_empty());
    paths.layered_dusty_gas = links;
    native::step_candidate_with_paths(parcels, paths, duration, endpoints, controls, cancelled)
}

const T: f64 = 400.0;
const P: f64 = 100_000.0;
const AREA: f64 = 0.01;

fn medium(permeability: f64) -> LayeredDustyGasMedium {
    let layers = [0.003, 0.007]
        .into_iter()
        .map(|length| {
            PorousGasLayer::try_new(
                IntrinsicPermeability::isotropic(permeability).unwrap(),
                length,
                0.0,
                "synthetic layered permeability",
            )
            .unwrap()
        })
        .collect();
    let barrier = DarcyGasBarrier::try_new(
        layers,
        [1.0, 0.0, 0.0],
        PorousGasSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (1000.0, 1e7),
            maximum_superficial_speed_m_s: 100.0,
        },
    )
    .unwrap();
    let pores = [(5e-8, 0.25, 1.5), (2e-7, 0.55, 3.0)]
        .into_iter()
        .map(|(radius, porosity, tortuosity)| {
            DustyGasPores::try_new(
                radius,
                porosity,
                tortuosity,
                "synthetic distinct pore layers",
            )
            .unwrap()
        })
        .collect();
    LayeredDustyGasMedium::try_new(barrier, pores).unwrap()
}
fn energy() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn solve() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.05,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    }
}
fn endpoints(gases: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    gases.keys().map(|&id| (id, energy())).collect()
}
fn pure(pressure: f64) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species("A", 0.028, 3.5, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025)],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(model, point(T, pressure), 0.001).unwrap()
}
fn binary(x: f64, pressure: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0 - x], 0.0), point(T, pressure), 0.001).unwrap()
}
fn ternary(x: Vec<f64>) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.028, 3.5, 0.0, 2e6),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, -1e6),
                species_with_energy_reference("C", 0.044, 4.0, 0.0, 4e6),
            ],
            Composition::new(CompositionBasis::MoleFraction, x).unwrap(),
        )
        .unwrap(),
        vec![
            transport("A", 1.7e-5, 0.025),
            transport("B", 2e-5, 0.15),
            transport("C", 1.5e-5, 0.02),
        ],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(model, point(T, P), 0.001).unwrap()
}
fn card(value: f64, pressures: (f64, f64)) -> MaterialCard {
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: BINARY_DIFFUSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    typed(QuantityKind::AbsoluteTemperature),
                    200.0,
                    1800.0,
                )
                .with_quantity(
                    "pressure",
                    typed(QuantityKind::Pressure),
                    pressures.0,
                    pressures.1,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic laminate diffusivity".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gas pair".into(),
            phase: "gas".into(),
            process: "finite laminate regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn binding<'a>(
    names: [&str; 2],
    card: &'a MaterialCard,
    context: &'a QueryPoint,
) -> BinaryDiffusionMaterial<'a> {
    BinaryDiffusionMaterial {
        species: names.map(|name| SpeciesId::new(name).unwrap()),
        material: card,
        context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn link<'a>(
    medium: &'a LayeredDustyGasMedium,
    pairs: &'a [BinaryDiffusionMaterial<'a>],
    context: &'a QueryPoint,
) -> LayeredDustyGasNetworkLink<'a> {
    LayeredDustyGasNetworkLink {
        id: 7,
        parcels: [10, 20],
        medium,
        area_m2: AREA,
        pairs,
        gas_context: context,
        solve: DustyGasFaceSolve {
            maximum_temperature_difference_k: 5.0,
            maximum_relative_pressure_difference: 0.5,
            gradient_tolerance_per_m: 1e-8,
        },
    }
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:e} != {expected:e}; tolerance {absolute:e}"
    );
}

#[test]
fn heterogeneous_knudsen_resistance_changes_actual_inventory_and_carried_energy() {
    let gases = BTreeMap::from([(10, pure(110_000.0)), (20, pure(90_000.0))]);
    let saved = gases.clone();
    let bounds = endpoints(&gases);
    let context = point(T, P);
    let stack = medium(0.0);
    let edge = link(&stack, &[], &context);
    let dt = 0.01;
    let candidate = step_candidate(&gases, &[edge], dt, &bounds, solve(), || false).unwrap();
    // Independent pure-gas serial wall resistance, retaining each pore geometry.
    let speed = (8.0 * R * T / (core::f64::consts::PI * 0.028)).sqrt();
    let dk = [
        (0.25 / 1.5) * (2.0 * 5e-8 / 3.0) * speed,
        (0.55 / 3.0) * (2.0 * 2e-7 / 3.0) * speed,
    ];
    let resistance = 0.003 / dk[0] + 0.007 / dk[1];
    let expected_mass = dt * AREA * 0.028 * (110_000.0 - 90_000.0) / (R * T * resistance);
    near(
        gases[&10].mass_kg() - candidate.parcels()[&10].mass_kg(),
        expected_mass,
        3e-17,
    );
    near(
        candidate.parcels()[&20].mass_kg() - gases[&20].mass_kg(),
        expected_mass,
        3e-17,
    );
    let native = candidate.step();
    let stream = &native.exchange.fluxes[0];
    let enthalpy = 3.5 * R * T / 0.028;
    near(stream.source_enthalpy_j_kg, enthalpy, enthalpy * 3e-12);
    near(
        stream.carried_enthalpy_j,
        expected_mass * enthalpy,
        expected_mass * enthalpy * 3e-11,
    );
    assert_eq!(native.exchange.fluxes.len(), 1);
    assert_eq!(
        candidate.step().layered_dusty_gas[0].flux.species[0].viscous_molar_flux_mol_m2_s,
        0.0
    );
    assert!(native.exchange.total_energy_residual_j.abs() <= solve().exchange.energy_tolerance_j);
    assert!(native.entropy_change_j_k > 0.0);
    assert!(candidate.parcels()[&10].state().state().pressure < 110_000.0);
    assert!(candidate.parcels()[&20].state().state().pressure > 90_000.0);
    for id in gases.keys() {
        assert_eq!(candidate.parcels()[id].volume_m3(), gases[id].volume_m3());
    }
    // The original Darcy contribution joins the same inventory exactly once.
    let permeable = medium(2e-14);
    let with_viscous = step_candidate(
        &gases,
        &[link(&permeable, &[], &context)],
        dt,
        &bounds,
        solve(),
        || false,
    )
    .unwrap();
    assert!(
        with_viscous.step().layered_dusty_gas[0].flux.species[0].viscous_molar_flux_mol_m2_s > 0.0
    );
    assert!(with_viscous.parcels()[&10].mass_kg() < candidate.parcels()[&10].mass_kg());
    assert_eq!(gases, saved);
}

#[test]
fn multicomponent_laminate_reuses_original_flux_and_reference_energy_update() {
    // fs-thermochem mixture.rs:349 sorts the fractions into canonical SpeciesId order (A, B, C)
    // and :366-367 refuse unless their `Iterator::sum` is bit-exact 1.0 (d571a4691d);
    // 0.70 + 0.20 + 0.10 sums to 1 - 2^-53. C closes the sum as 1 - s, with s = A + B taken
    // the same way. With the partial sum s in [0.5, 2], Sterbenz's lemma makes 1 - s exact,
    // so s + (1 - s) == 1 exactly. Here s = 0.8999999999999999 and C = 0.10000000000000009.
    let leading: [f64; 2] = [0.70, 0.20];
    let others: f64 = leading.iter().sum();
    let fractions = vec![leading[0], leading[1], 1.0 - others];
    assert_eq!(fractions.iter().sum::<f64>().to_bits(), 1.0_f64.to_bits());
    let gases = BTreeMap::from([
        (10, ternary(fractions)),
        (20, ternary(vec![0.15, 0.25, 0.60])),
        (30, ternary(vec![0.15, 0.25, 0.60])),
    ]);
    let saved = gases.clone();
    let bounds = endpoints(&gases);
    let context = point(T, P);
    let cards = [
        card(2e-5, (1000.0, 1e7)),
        card(1.5e-5, (1000.0, 1e7)),
        card(3e-5, (1000.0, 1e7)),
    ];
    let pairs = [
        binding(["A", "B"], &cards[0], &context),
        binding(["A", "C"], &cards[1], &context),
        binding(["B", "C"], &cards[2], &context),
    ];
    let stack = medium(2e-14);
    let edge = link(&stack, &pairs, &context);
    let dt = 0.01;
    let flux = gases[&10]
        .layered_dusty_gas_flux_to(
            &gases[&20],
            &stack,
            AREA,
            &pairs,
            &context,
            edge.solve,
            || false,
        )
        .unwrap();
    let transfers: Vec<_> = flux
        .species
        .iter()
        .filter_map(|row| {
            let mass = dt * row.mass_rate_left_to_right_kg_s;
            (mass != 0.0).then(|| SpeciesNetworkTransfer {
                donor: if mass > 0.0 { 10 } else { 20 },
                receiver: if mass > 0.0 { 20 } else { 10 },
                species: row.species.clone(),
                mass_kg: mass.abs(),
            })
        })
        .collect();
    let expected =
        network::exchange_candidate(&gases, &transfers, &bounds, solve().exchange, || false)
            .unwrap();
    let candidate = step_candidate(&gases, &[edge], dt, &bounds, solve(), || false).unwrap();
    assert_eq!(candidate.parcels(), expected.parcels());
    assert_eq!(&candidate.step().exchange, expected.step());
    assert_eq!(candidate.step().layered_dusty_gas[0].flux, flux);
    assert_eq!(candidate.step().layered_dusty_gas[0].flux.pairs.len(), 3);
    assert!(flux.species.iter().any(|row| row.molar_flux_mol_m2_s > 0.0));
    assert!(flux.species.iter().any(|row| row.molar_flux_mol_m2_s < 0.0));
    assert!(flux.net_molar_flux_mol_m2_s.abs() > 1e-12);
    let native = candidate.step();
    assert!(native.exchange.maximum_species_mass_residual_kg <= solve().exchange.mass_tolerance_kg);
    assert!(native.exchange.total_energy_residual_j.abs() <= solve().exchange.energy_tolerance_j);
    assert!(native.entropy_change_j_k > 0.0);
    assert_eq!(candidate.parcels()[&30], gases[&30]);
    assert_eq!(gases, saved);
}

#[test]
fn every_laminate_reads_one_snapshot_and_later_steps_transport_arrivals() {
    let gases = BTreeMap::from([
        (10, binary(0.75, P)),
        (20, binary(0.25, P)),
        (30, binary(0.25, P)),
    ]);
    let bounds = endpoints(&gases);
    let stack = medium(0.0);
    let context = point(T, P);
    // Onward transport of an arrival needs molecular friction weak against B's
    // Knudsen back-flow. Step 1 is isobaric Graham flux, g = |N_B/N_A| = sqrt 7:
    // parcel 20 gains A, loses g - 1 = 1.65 times as many net moles, and cools
    // about 3 mK, tau/eps = -(g - 1)/(sum x_i cv_i/R) = -0.94 (eps the relative
    // A gain, tau the relative temperature change). In the later step B returns
    // 30 -> 20 and drags A with it, uphill (contracts/part-03.md:3210-3211),
    // unless D_AB exceeds
    //     D* = D0 R_AB x_A ((g - 1) - tau/eps) / (K_B (1 + x_A tau/eps)) = 8.36e-5 m2/s,
    // with K_B and R_AB the laminate's Knudsen and molecular resistances at any
    // D0 (D0 R_AB does not depend on D0); D* is 4.07e-5 if parcel 20 stayed
    // isothermal. An independent binary dusty-gas oracle gives link 8 in the
    // later step N_A = -1.01e-7 and N_B = -7.66e-7 mol/m2/s at 2e-5, and
    // N_A = +1.75e-7 and N_B = -3.86e-6 mol/m2/s at 2e-4, where 4.9e-13 kg of A
    // moves 20 -> 30. 2e-4 is a fixture input chosen 2.39 times above D*, not a
    // pass threshold and not a measured card: it sits at the high end of binary
    // gas diffusivities at 1 bar (light-gas pairs are about 0.6e-4 to 0.7e-4
    // m2/s at room temperature).
    let card = card(2e-4, (1000.0, 1e7));
    let pairs = [binding(["A", "B"], &card, &context)];
    let first = link(&stack, &pairs, &context);
    let second = LayeredDustyGasNetworkLink {
        id: 8,
        parcels: [20, 30],
        ..first
    };
    let candidate =
        step_candidate(&gases, &[first, second], 0.01, &bounds, solve(), || false).unwrap();
    assert!(
        candidate.step().layered_dusty_gas[1]
            .flux
            .species
            .iter()
            .all(|row| row.mass_rate_left_to_right_kg_s == 0.0)
    );
    assert_eq!(candidate.parcels()[&30], gases[&30]);
    let reordered = step_candidate(
        &gases,
        &[
            LayeredDustyGasNetworkLink { id: 2, ..second },
            LayeredDustyGasNetworkLink { id: 900, ..first },
        ],
        0.01,
        &bounds,
        solve(),
        || false,
    )
    .unwrap();
    assert_eq!(candidate.parcels(), reordered.parcels());
    assert_eq!(candidate.step().exchange, reordered.step().exchange);
    let next = step_candidate(
        candidate.parcels(),
        &[first, second],
        0.01,
        &bounds,
        solve(),
        || false,
    )
    .unwrap();
    let a = SpeciesId::new("A").unwrap();
    assert!(
        next.step()
            .exchange
            .fluxes
            .iter()
            .any(|row| row.transfer.donor == 20
                && row.transfer.receiver == 30
                && row.transfer.species == a
                && row.transfer.mass_kg > 0.0)
    );
    assert!(
        inventory(next.parcels()[&30].model(), next.parcels()[&30].mass_kg()).unwrap()[&a]
            > inventory(gases[&30].model(), gases[&30].mass_kg()).unwrap()[&a]
    );
}

#[test]
fn layered_and_existing_paths_share_one_gross_species_allowance() {
    let gases = BTreeMap::from([(10, pure(110_000.0)), (20, pure(90_000.0))]);
    let saved = gases.clone();
    let bounds = endpoints(&gases);
    let context = point(T, P);
    let stack = medium(0.0);
    let mut layered = link(&stack, &[], &context);
    let dt = 0.01;
    let homogeneous = DustyGasMedium::try_new(
        DarcyGasBarrier::try_new(
            vec![stack.barrier().layers()[0].clone()],
            [1.0, 0.0, 0.0],
            stack.barrier().support(),
        )
        .unwrap(),
        stack.pore_layers()[0].clone(),
    )
    .unwrap();
    let mut existing = DustyGasNetworkLink {
        id: 19,
        parcels: [10, 20],
        medium: &homogeneous,
        area_m2: AREA,
        pairs: &[],
        gas_context: &context,
        solve: layered.solve,
    };
    let layered_rate = gases[&10]
        .layered_dusty_gas_flux_to(
            &gases[&20],
            &stack,
            AREA,
            &[],
            &context,
            layered.solve,
            || false,
        )
        .unwrap()
        .species[0]
        .mass_rate_left_to_right_kg_s;
    let existing_rate = gases[&10]
        .dusty_gas_flux_to(
            &gases[&20],
            &homogeneous,
            AREA,
            &[],
            &context,
            existing.solve,
            || false,
        )
        .unwrap()
        .species[0]
        .mass_rate_left_to_right_kg_s;
    let target = 0.03 * gases[&10].mass_kg();
    layered.area_m2 *= target / (layered_rate * dt);
    existing.area_m2 *= target / (existing_rate * dt);
    step_candidate(&gases, &[layered], dt, &bounds, solve(), || false).unwrap();
    native::step_candidate_with_paths(
        &gases,
        GasTransportPaths {
            dusty_gas: &[existing],
            ..GasTransportPaths::default()
        },
        dt,
        &bounds,
        solve(),
        || false,
    )
    .unwrap();
    assert!(matches!(
        step_candidate_with_paths(
            &gases,
            GasTransportPaths {
                dusty_gas: &[existing],
                ..GasTransportPaths::default()
            },
            &[layered],
            dt,
            &bounds,
            solve(),
            || false
        ),
        Err(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { id: 10, .. }
        ))
    ));
    assert_eq!(gases, saved);
}

#[test]
fn inactive_faces_sources_and_late_cancellation_preserve_original_gases() {
    let gases = BTreeMap::from([(10, binary(0.75, P)), (20, binary(0.25, P))]);
    let saved = gases.clone();
    let bounds = endpoints(&gases);
    let context = point(T, P);
    let stack = medium(0.0);
    let source = card(2e-5, (1000.0, 1e7));
    let pairs = [binding(["A", "B"], &source, &context)];
    let edge = link(&stack, &pairs, &context);
    let closed = LayeredDustyGasNetworkLink {
        area_m2: 0.0,
        pairs: &[],
        ..edge
    };
    let idle = step_candidate(&gases, &[closed], 0.01, &bounds, solve(), || false).unwrap();
    assert_eq!(idle.parcels(), &gases);
    assert!(idle.step().layered_dusty_gas[0].flux.face_state.is_none());
    assert!(idle.step().exchange.fluxes.is_empty());
    let empty = step_candidate(&gases, &[], 0.01, &bounds, solve(), || false).unwrap();
    assert_eq!(empty.parcels(), &gases);
    assert!(matches!(
        step_candidate(
            &gases,
            &[LayeredDustyGasNetworkLink { pairs: &[], ..edge }],
            0.01,
            &bounds,
            solve(),
            || false
        ),
        Err(GasTransportError::DustyGas(DustyGasNetworkError::Face {
            id: 7,
            ..
        }))
    ));
    // Mean pressure is covered; actual original endpoint pressures are not.
    let narrow = card(2e-5, (99_000.0, 101_000.0));
    let narrow_pairs = [binding(["A", "B"], &narrow, &context)];
    let outside = BTreeMap::from([(10, binary(0.75, 110_000.0)), (20, binary(0.25, 90_000.0))]);
    assert!(matches!(
        step_candidate(
            &outside,
            &[link(&stack, &narrow_pairs, &context)],
            0.01,
            &endpoints(&outside),
            solve(),
            || false
        ),
        Err(GasTransportError::DustyGas(DustyGasNetworkError::Face {
            id: 7,
            ..
        }))
    ));
    let mut polls = 0;
    step_candidate(&gases, &[edge], 0.01, &bounds, solve(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut calls = 0;
        assert_eq!(
            step_candidate(&gases, &[edge], 0.01, &bounds, solve(), || {
                calls += 1;
                calls == stop
            }),
            Err(GasTransportError::Cancelled)
        );
    }
    assert_eq!(gases, saved);
}

fn matched_binary(x: f64) -> ClosedGasParcel {
    // Equal molar masses and Cp close an independent linear composition mode.
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species("A", 0.028, 3.5, 0.0), species("B", 0.028, 3.5, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(model, point(T, P), 0.001).unwrap()
}
fn adaptive(tolerance: f64) -> EffusionAdvanceControls {
    EffusionAdvanceControls {
        initial_step_s: 200.0,
        minimum_step_s: 1e-8,
        maximum_step_s: 200.0,
        maximum_attempts: 20_000,
        maximum_accepted_intervals: 4096,
        relative_tolerance: 0.0,
        temperature_tolerance_k: 1e-6,
        pressure_tolerance_pa: 1e-3,
        species_mass_tolerance_kg: tolerance,
        stream_energy_tolerance_j: 1e-6,
    }
}

#[test]
fn closed_laminate_and_molecular_paths_refine_to_independent_binary_relaxation() {
    let gases = BTreeMap::from([(10, matched_binary(0.75)), (20, matched_binary(0.25))]);
    let saved = gases.clone();
    let bounds = endpoints(&gases);
    let context = point(T, P);
    let stack = medium(0.0);
    let source = card(2e-5, (1000.0, 1e7));
    let pairs = [binding(["A", "B"], &source, &context)];
    let layers = [link(&stack, &pairs, &context)];
    // This molecular face is a separate exposed path, not a second model of
    // the same laminate aperture.
    let molecular = [DiffusionNetworkLink {
        id: 19,
        parcels: [10, 20],
        pairs: &pairs,
        geometry: GasDiffusionGeometry {
            area_m2: 0.001,
            distance_m: 0.005,
        },
        solve: DiffusionFaceSolve {
            maximum_temperature_difference_k: 0.01,
            maximum_pressure_difference_pa: 1e-3,
            gradient_tolerance_per_m: 1e-8,
            molar_flux_tolerance_mol_m2_s: 1e-10,
        },
    }];
    let paths = GasTransportPaths {
        layered_dusty_gas: &layers,
        diffusion: &molecular,
        ..Default::default()
    };
    let speed = (8.0 * R * T / (core::f64::consts::PI * 0.028)).sqrt();
    let rk = 0.003 / ((0.25 / 1.5) * (2.0 * 5e-8 / 3.0) * speed)
        + 0.007 / ((0.55 / 3.0) * (2.0 * 2e-7 / 3.0) * speed);
    let rm = 0.003 / ((0.25 / 1.5) * 2e-5) + 0.007 / ((0.55 / 3.0) * 2e-5);
    let rate = (AREA / (rk + rm) + 0.001 * 2e-5 / 0.005) / 0.001;
    let exact_gap = 0.5 * (-2.0 * rate * 200.0).exp();
    let mut errors = Vec::new();
    let mut counts = Vec::new();
    for tolerance in [2e-7, 2e-8, 2e-9] {
        let result = advance_candidate(
            &gases,
            paths,
            &[],
            &[],
            &[],
            200.0,
            &bounds,
            solve(),
            adaptive(tolerance),
            || false,
        )
        .unwrap();
        let x = |id| {
            result.parcels()[&id]
                .model()
                .thermodynamics()
                .mole_composition()
                .fractions()[0]
        };
        errors.push(((x(10) - x(20)) - exact_gap).abs());
        counts.push(result.step().accepted_intervals);
        assert!(result.step().rejected_intervals > 0);
        assert!(result.step().maximum_accepted_error <= 1.0);
        assert_eq!(
            result.step().accepted_steps.len(),
            2 * result.step().accepted_intervals
        );
        assert!(result.step().boundaries.is_empty());
        assert!(
            result.step().balance.energy_residual_j.abs() <= solve().exchange.energy_tolerance_j
        );
        assert!(
            result.step().balance.accumulated_energy_residual_j.abs()
                <= solve().exchange.energy_tolerance_j
        );
        let audit =
            native::audit_closed_window(&gases, result.parcels(), solve(), || false).unwrap();
        assert!(audit.entropy_change_j_k > 0.0);
        assert!(audit.maximum_species_mass_residual_kg <= solve().exchange.mass_tolerance_kg);
        near(
            result
                .step()
                .accepted_steps
                .iter()
                .map(|s| s.duration_s)
                .sum(),
            200.0,
            1e-9,
        );
        let a = SpeciesId::new("A").unwrap();
        let mut transferred_a = 0.0;
        for step in &result.step().accepted_steps {
            assert_eq!(step.transport.layered_dusty_gas.len(), 1);
            assert_eq!(step.transport.diffusion.len(), 1);
            for stream in &step.exchange.fluxes {
                if stream.transfer.species == a {
                    transferred_a += if stream.transfer.donor == 10 {
                        stream.transfer.mass_kg
                    } else {
                        -stream.transfer.mass_kg
                    };
                }
            }
        }
        near(
            gases[&10].species_mass_kg(&a).unwrap()
                - result.parcels()[&10].species_mass_kg(&a).unwrap(),
            transferred_a,
            solve().exchange.mass_tolerance_kg,
        );
        for gas in result.parcels().values() {
            near(gas.state().state().temperature, T, 1e-6);
            near(gas.state().state().pressure, P, 1e-3);
        }
    }
    assert!(errors[2] < errors[1] && errors[1] < errors[0], "{errors:?}");
    assert!(errors[2] < 1e-3, "{errors:?}");
    assert!(counts[2] > counts[1] && counts[1] > counts[0], "{counts:?}");
    assert_eq!(gases, saved);
}

#[test]
fn adaptive_laminate_work_accuracy_source_and_late_cancellation_are_atomic() {
    let gases = BTreeMap::from([(10, matched_binary(0.75)), (20, matched_binary(0.25))]);
    let saved = gases.clone();
    let bounds = endpoints(&gases);
    let context = point(T, P);
    let stack = medium(0.0);
    let source = card(2e-5, (1000.0, 1e7));
    let pairs = [binding(["A", "B"], &source, &context)];
    let layers = [link(&stack, &pairs, &context)];
    let paths = GasTransportPaths {
        layered_dusty_gas: &layers,
        ..Default::default()
    };
    assert!(matches!(
        advance_candidate(
            &gases,
            paths,
            &[],
            &[],
            &[],
            200.0,
            &bounds,
            solve(),
            EffusionAdvanceControls {
                maximum_attempts: 3,
                ..adaptive(2e-9)
            },
            || false
        ),
        Err(EffusionAdvanceError::WorkBudget {
            attempted_steps: 3,
            ..
        })
    ));
    assert!(matches!(
        advance_candidate(
            &gases,
            paths,
            &[],
            &[],
            &[],
            1.0,
            &bounds,
            solve(),
            EffusionAdvanceControls {
                initial_step_s: 1.0,
                minimum_step_s: 1.0,
                maximum_step_s: 1.0,
                ..adaptive(1e-16)
            },
            || false
        ),
        Err(EffusionAdvanceError::MinimumStep { .. })
    ));
    let unsupported = card(2e-5, (110_000.0, 120_000.0));
    let bad_pairs = [binding(["A", "B"], &unsupported, &context)];
    let bad_layers = [link(&stack, &bad_pairs, &context)];
    assert!(matches!(
        advance_candidate(
            &gases,
            GasTransportPaths {
                layered_dusty_gas: &bad_layers,
                ..Default::default()
            },
            &[],
            &[],
            &[],
            200.0,
            &bounds,
            solve(),
            EffusionAdvanceControls {
                maximum_attempts: 3,
                ..adaptive(2e-9)
            },
            || false
        ),
        Err(EffusionAdvanceError::Step {
            error: GasTransportError::DustyGas(DustyGasNetworkError::Face { id: 7, .. }),
            ..
        })
    ));
    let mut polls = 0;
    advance_candidate(
        &gases,
        paths,
        &[],
        &[],
        &[],
        20.0,
        &bounds,
        solve(),
        adaptive(2e-7),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut calls = 0;
    assert_eq!(
        advance_candidate(
            &gases,
            paths,
            &[],
            &[],
            &[],
            20.0,
            &bounds,
            solve(),
            adaptive(2e-7),
            || {
                calls += 1;
                calls == polls
            }
        ),
        Err(EffusionAdvanceError::Cancelled)
    );
    assert_eq!(gases, saved);
}

#[test]
fn inactive_layered_faces_count_against_the_original_adaptive_history_bound() {
    let gases = BTreeMap::from([(10, matched_binary(0.5)), (20, matched_binary(0.5))]);
    let bounds = endpoints(&gases);
    let context = point(T, P);
    let stack = medium(0.0);
    let layers = [LayeredDustyGasNetworkLink {
        area_m2: 0.0,
        ..link(&stack, &[], &context)
    }];
    let paths = GasTransportPaths {
        layered_dusty_gas: &layers,
        ..Default::default()
    };
    // Each accepted report retains 1 step + 2 parcel accounts + 1 inactive
    // laminate. Both half steps must fit; the former count omitted the laminate.
    let limited = advance_candidate_with_row_limit(
        &gases,
        paths,
        &[],
        &[],
        &[],
        1.0,
        &bounds,
        solve(),
        adaptive(1e-6),
        6,
        || false,
    );
    assert!(matches!(
        limited,
        Err(EffusionAdvanceError::WorkBudget {
            attempted_steps: 3,
            accepted_intervals: 0,
            ..
        })
    ));
    let complete = advance_candidate_with_row_limit(
        &gases,
        paths,
        &[],
        &[],
        &[],
        1.0,
        &bounds,
        solve(),
        adaptive(1e-6),
        8,
        || false,
    )
    .unwrap();
    assert_eq!(complete.parcels(), &gases);
    assert_eq!(complete.step().accepted_steps.len(), 2);
    assert!(
        complete
            .step()
            .accepted_steps
            .iter()
            .all(|s| s.transport.layered_dusty_gas.len() == 1)
    );
}
