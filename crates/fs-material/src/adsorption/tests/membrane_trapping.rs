//! Source-backed mobile/trap coupling and vacuum-domain admission.
use super::*;
#[path = "membrane_trapping/adaptive.rs"]
mod adaptive;
use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::{
    MEMBRANE_DIFFUSIVITY_DIMS, MOLAR_HENRY_SOLUBILITY_DIMS, MembranePermeabilityMaterial,
    transient::{
        MEMBRANE_DIFFUSIVITY_PROPERTY, MOLAR_HENRY_SOLUBILITY_PROPERTY,
        MembraneDiffusionLayerSource, ResolvedMembraneDiffusion,
        slab::trapping::TrappingMembraneSolve,
        trapping::{
            LangmuirMembraneTrapSource, ResolvedTrappingMembraneDiffusion,
            ResolvedTrappingMembraneError,
        },
    },
};

fn pure_catalog() -> ThermochemicalGasModel {
    // The active-species NASA owner rejects zero fractions. The original gas
    // composition adapter owns removal of exactly inactive catalog species.
    model(vec![0.5, 0.5], 0.0)
        .with_composition(
            fs_thermochem::Composition::new(
                fs_thermochem::CompositionBasis::MoleFraction,
                vec![1.0, 0.0],
            )
            .unwrap(),
        )
        .unwrap()
}

fn source_claim(name: &str, dims: Dims, value: f64, upper: f64) -> PropertyClaim {
    let mut property = claim(name, dims, value);
    property.validity = ValidityDomain::unconstrained()
        .with_quantity("temperature", temp_quantity(), 200.0, 800.0)
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            0.0,
            upper,
        );
    property
}
fn vacuum_source(upper: f64) -> AdsorptionSource {
    AdsorptionSource {
        material: card(
            "vacuum-qualified trap constituent",
            vec![
                source_claim(
                    LANGMUIR_CAPACITY_PROPERTY,
                    LANGMUIR_CAPACITY_DIMS,
                    0.5,
                    upper,
                ),
                source_claim(
                    LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
                    LANGMUIR_AFFINITY_DIMS,
                    2e-5,
                    upper,
                ),
                source_claim(
                    LANGMUIR_DESORPTION_RATE_PROPERTY,
                    LANGMUIR_RATE_DIMS,
                    0.2,
                    upper,
                ),
            ],
        ),
        point: point(T0, PRESSURE),
        selection: MaterialPropertySelection::SingleClaimOnly,
        reference_temperature_k: T0,
    }
}
fn law(catalog: &ThermochemicalGasModel, source: AdsorptionSource) -> Arc<LangmuirAdsorbent> {
    Arc::new(
        LangmuirAdsorbent::try_new(
            catalog,
            SpeciesId::new("A").unwrap(),
            source,
            calorics(catalog, 0.0, HEAT),
        )
        .unwrap(),
    )
}
fn membrane_card() -> MaterialCard {
    card(
        "mobile polymer constituent",
        vec![
            source_claim(
                MEMBRANE_DIFFUSIVITY_PROPERTY,
                MEMBRANE_DIFFUSIVITY_DIMS,
                1e-10,
                1e7,
            ),
            source_claim(
                MOLAR_HENRY_SOLUBILITY_PROPERTY,
                MOLAR_HENRY_SOLUBILITY_DIMS,
                2e-5,
                1e7,
            ),
        ],
    )
}
fn mobile<'a>(card: &'a MaterialCard, context: &'a QueryPoint) -> ResolvedMembraneDiffusion<'a> {
    let layer = MembraneDiffusionLayerSource {
        thickness_m: 1e-4,
        cells: 8,
        material: MembranePermeabilityMaterial {
            species: SpeciesId::new("A").unwrap(),
            material: card,
            context,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    };
    ResolvedMembraneDiffusion::new(T0, 0.01, &[layer], &[], &[0.0; 8], [PRESSURE, 0.0], || {
        false
    })
    .unwrap()
}

#[test]
fn vacuum_affinity_is_source_qualified_without_admitting_a_vacuum_gas_state() {
    let catalog = pure_catalog();
    let original = adsorbent(&catalog, 2e-5, 0.2);
    assert!(matches!(
        original.properties(T0, 0.0),
        Err(AdsorptionError::Material(_))
    ));
    let qualified = law(&catalog, vacuum_source(1e7));
    let vacuum = qualified.properties(T0, 0.0).unwrap();
    let positive = qualified.properties(T0, PRESSURE).unwrap();
    assert_eq!(vacuum.affinity_pa_inverse, positive.affinity_pa_inverse);
    assert_eq!(
        vacuum.desorption_rate_s_inverse,
        positive.desorption_rate_s_inverse
    );
    assert_eq!(vacuum.source.state_coordinate("pressure"), Some(0.0));
    assert!(qualified.properties(T0, -1.0).is_err());
    assert!(ClosedGasParcel::new(catalog, point(T0, 0.0), 0.001).is_err());
}

#[test]
fn original_trap_and_mobile_sources_share_the_same_conservative_native_step() {
    let catalog = pure_catalog();
    let trap = law(&catalog, vacuum_source(PRESSURE));
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let state = ResolvedTrappingMembraneDiffusion::new(
        mobile(&card, &context),
        &[Some(LangmuirMembraneTrapSource {
            adsorbent: trap.clone(),
            dry_skeleton_density_kg_m3: 4.0,
        })],
        &[0.0; 8],
        [PRESSURE, 0.0],
        || false,
    )
    .unwrap();
    let mut work = state.workspace().unwrap();
    let controls = TrappingMembraneSolve::default();
    let accepted = state
        .step_candidate([PRESSURE, 0.0], 10.0, controls, &mut work, || false)
        .unwrap();
    let mut raw_work = state.state().workspace().unwrap();
    let raw = state
        .state()
        .step_candidate([PRESSURE, 0.0], 10.0, controls, &mut raw_work, || false)
        .unwrap();
    assert_eq!(accepted.state().state(), raw.state());
    assert_eq!(accepted.step(), raw.step());
    assert!(accepted.step().trapped_change_mol > 0.0);
    for support in accepted.state().trap_source_support()[0].as_ref().unwrap() {
        assert_eq!(
            support.source.card_identity(),
            trap.source.material.content_hash()
        );
        assert_eq!(support.affinity_pa_inverse, 2e-5);
        assert_eq!(support.desorption_rate_s_inverse, 0.2);
    }
    let saved = state.state().clone();
    let receipts = state.trap_source_support().to_vec();
    // The polymer supports this pressure; the original trap constituent does not.
    assert!(matches!(
        state.step_candidate([2.0 * PRESSURE, 0.0], 10.0, controls, &mut work, || false),
        Err(ResolvedTrappingMembraneError::Adsorption(
            AdsorptionError::Material(_)
        ))
    ));
    assert!(matches!(
        state.step_candidate([PRESSURE, 0.0], 10.0, controls, &mut work, || true),
        Err(ResolvedTrappingMembraneError::Cancelled)
    ));
    assert_eq!(state.state(), &saved);
    assert_eq!(state.trap_source_support(), receipts);
}

#[test]
fn different_humidity_or_unowned_mixture_coordinates_refuse_composite_binding() {
    use fs_qty::semantic::{CompositionBasis, MoistureBasis};
    let catalog = pure_catalog();
    let card = membrane_card();
    let humidity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::RelativeHumidity),
        ValueForm::Static,
    ));
    let mut source = vacuum_source(PRESSURE);
    source.point = source
        .point
        .with_quantity("humidity", humidity, 0.5)
        .unwrap();
    let trap = law(&catalog, source);
    let context = point(T0, PRESSURE)
        .with_quantity("humidity", humidity, 0.4)
        .unwrap();
    assert!(matches!(
        ResolvedTrappingMembraneDiffusion::new(
            mobile(&card, &context),
            &[Some(LangmuirMembraneTrapSource {
                adsorbent: trap,
                dry_skeleton_density_kg_m3: 4.0
            })],
            &[0.0; 8],
            [PRESSURE, 0.0],
            || false
        ),
        Err(ResolvedTrappingMembraneError::Invalid(_))
    ));

    let fraction = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Composition(CompositionBasis::MoleFraction),
        ValueForm::Static,
    ));
    let mut source = vacuum_source(PRESSURE);
    source.point = source
        .point
        .with_quantity("mixture_fraction", fraction, 0.5)
        .unwrap();
    let trap = law(&catalog, source);
    let context = point(T0, PRESSURE)
        .with_quantity("mixture_fraction", fraction, 0.5)
        .unwrap();
    assert!(matches!(
        ResolvedTrappingMembraneDiffusion::new(
            mobile(&card, &context),
            &[Some(LangmuirMembraneTrapSource {
                adsorbent: trap,
                dry_skeleton_density_kg_m3: 4.0
            })],
            &[0.0; 8],
            [PRESSURE, 0.0],
            || false
        ),
        Err(ResolvedTrappingMembraneError::Invalid(_))
    ));
}

#[test]
fn finite_binding_qualifies_actual_baths_and_preserves_all_native_inventories() {
    use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::{
        slab::finite::MolarPressureCapacity, trapping::ResolvedFiniteTrappingMembraneDiffusion,
    };
    let catalog = pure_catalog();
    let trap = law(&catalog, vacuum_source(PRESSURE));
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let film = ResolvedTrappingMembraneDiffusion::new(
        mobile(&card, &context),
        &[Some(LangmuirMembraneTrapSource {
            adsorbent: trap,
            dry_skeleton_density_kg_m3: 4.0,
        })],
        &[0.0; 8],
        [PRESSURE, 0.0],
        || false,
    )
    .unwrap();
    let capacity = MolarPressureCapacity::try_new(2e-11).unwrap();
    assert!(matches!(
        ResolvedFiniteTrappingMembraneDiffusion::new(
            film.clone(),
            [capacity; 2],
            [2.0 * PRESSURE, 0.0],
            || false
        ),
        Err(ResolvedTrappingMembraneError::Adsorption(
            AdsorptionError::Material(_)
        ))
    ));
    let state = ResolvedFiniteTrappingMembraneDiffusion::new(
        film,
        [capacity; 2],
        [0.5 * PRESSURE, 0.0],
        || false,
    )
    .unwrap();
    assert_eq!(
        state.source_support().layers[0][1].state_coordinate("pressure"),
        Some(0.5 * PRESSURE)
    );
    assert_eq!(
        state.trap_source_support()[0].as_ref().unwrap()[1]
            .source
            .state_coordinate("pressure"),
        Some(0.5 * PRESSURE)
    );
    let total = state.state().total_amount_mol().unwrap();
    let controls = TrappingMembraneSolve::default();
    let mut work = state.workspace().unwrap();
    let candidate = state
        .step_candidate(10.0, controls, &mut work, || false)
        .unwrap();
    let mut native_work = state.state().workspace().unwrap();
    let native = state
        .state()
        .step_candidate(10.0, controls, &mut native_work, || false)
        .unwrap();
    assert_eq!(candidate.state().state(), native.state());
    assert_eq!(candidate.step(), native.step());
    assert!((candidate.state().state().total_amount_mol().unwrap() - total).abs() < 1e-13);
    assert!(candidate.state().reservoir_pressures_pa()[0] < state.reservoir_pressures_pa()[0]);
    assert!(candidate.state().reservoir_pressures_pa()[1] > 0.0);
    assert!(candidate.step().trapped_change_mol > 0.0);
    // Accepted support includes the initial inlet, not only its depleted value.
    assert_eq!(
        candidate.state().source_support().layers[0][1].state_coordinate("pressure"),
        Some(0.5 * PRESSURE)
    );
}

#[test]
fn finite_trap_release_cannot_publish_an_unsupported_pressure_or_partial_inventory() {
    use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::{
        slab::finite::MolarPressureCapacity, trapping::ResolvedFiniteTrappingMembraneDiffusion,
    };
    let catalog = pure_catalog();
    let mut source = vacuum_source(3_500.0);
    source.point = point(T0, 0.0);
    let trap = law(&catalog, source);
    let card = membrane_card();
    let context = point(T0, 0.0);
    let layer = MembraneDiffusionLayerSource {
        thickness_m: 1e-4,
        cells: 1,
        material: MembranePermeabilityMaterial {
            species: SpeciesId::new("A").unwrap(),
            material: &card,
            context: &context,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    };
    let film = ResolvedMembraneDiffusion::new(T0, 0.01, &[layer], &[], &[0.0], [0.0; 2], || false)
        .unwrap();
    let trapped = ResolvedTrappingMembraneDiffusion::new(
        film,
        &[Some(LangmuirMembraneTrapSource {
            adsorbent: trap,
            dry_skeleton_density_kg_m3: 4.0,
        })],
        &[(0.01 * 1e-4) * 2.0],
        [0.0; 2],
        || false,
    )
    .unwrap();
    let capacity = MolarPressureCapacity::try_new(2e-11).unwrap();
    let state =
        ResolvedFiniteTrappingMembraneDiffusion::new(trapped, [capacity; 2], [0.0; 2], || false)
            .unwrap();
    let saved = state.state().clone();
    let support = state.source_support().clone();
    let traps = state.trap_source_support().to_vec();
    let mut work = state.workspace().unwrap();
    let controls = TrappingMembraneSolve::default();
    let native = state
        .state()
        .step_candidate(1_000.0, controls, &mut work, || false)
        .unwrap();
    assert!(native.state().membrane().mobile().partial_pressures_pa()[0] > 3_500.0);
    assert!(matches!(
        state.step_candidate(1_000.0, controls, &mut work, || false),
        Err(ResolvedTrappingMembraneError::Adsorption(
            AdsorptionError::Material(_)
        ))
    ));
    assert!(matches!(
        state.step_candidate(1.0, controls, &mut work, || true),
        Err(ResolvedTrappingMembraneError::Cancelled)
    ));
    assert_eq!(state.state(), &saved);
    assert_eq!(state.source_support(), &support);
    assert_eq!(state.trap_source_support(), traps);
}
