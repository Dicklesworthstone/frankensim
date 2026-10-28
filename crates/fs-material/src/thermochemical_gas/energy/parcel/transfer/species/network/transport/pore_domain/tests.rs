//! Actual ideal-mixture/DGM inventories, with a separate linear-network oracle.
use super::*;
use crate::porous::permeability::{
    DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
};
use crate::state_point::MaterialPropertySelection;
use crate::thermochemical_gas::{
    energy::parcel::transfer::species::diffusion::{
        BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, dusty_gas::DustyGasPores,
    },
    energy::parcel::transfer::species::network::SpeciesNetworkSolve,
    tests::{point, species, transport},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};

const T: f64 = 400.0;
const P: f64 = 100_000.0;
const MASS: f64 = 0.028;
const D: f64 = 2e-5;
const BODY_AREA: f64 = 0.1;
const THROAT_AREA: f64 = 0.01;

fn medium(radius: f64, porosity: f64, tortuosity: f64, length: f64) -> DustyGasMedium {
    let barrier = DarcyGasBarrier::try_new(
        vec![
            PorousGasLayer::try_new(
                IntrinsicPermeability::isotropic(0.0).unwrap(),
                length,
                0.0,
                "synthetic finite pore-body/throat fixture",
            )
            .unwrap(),
        ],
        [1.0, 0.0, 0.0],
        PorousGasSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (1000.0, 1e7),
            maximum_superficial_speed_m_s: 100.0,
        },
    )
    .unwrap();
    DustyGasMedium::try_new(
        barrier,
        DustyGasPores::try_new(
            radius,
            porosity,
            tortuosity,
            "synthetic pore geometry; no measured authority",
        )
        .unwrap(),
    )
    .unwrap()
}
fn gas_model(x: f64) -> ThermochemicalGasModel {
    // Two distinguishable gases with equal molecular mass and caloric laws:
    // equimolar counterdiffusion leaves each node's total mass, T and p fixed.
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species("A", MASS, 3.5, 0.0), species("B", MASS, 3.5, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 1.7e-5, 0.025)],
        None,
    )
    .unwrap()
}
fn gas(x: f64, volume: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(gas_model(x), point(T, P), volume).unwrap()
}
fn x(parcel: &ClosedGasParcel) -> f64 {
    let law = parcel.model().thermodynamics();
    let a = SpeciesId::new("A").unwrap();
    let i = law
        .components()
        .iter()
        .position(|c| c.species() == &a)
        .unwrap();
    law.mole_composition().fractions()[i]
}
fn a_mass(parcel: &ClosedGasParcel) -> f64 {
    parcel.mass_kg() * x(parcel)
}
fn coefficient() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: D,
                dims: BINARY_DIFFUSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                200.0,
                1800.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic common molecular diffusivity".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "A/B equal-mass fixture".into(),
            phase: "gas".into(),
            process: "pore-network analytic oracle".into(),
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
fn face<'a>(
    medium: &'a DustyGasMedium,
    area_m2: f64,
    pairs: &'a [BinaryDiffusionMaterial<'a>],
    context: &'a QueryPoint,
) -> PoreFace<'a> {
    PoreFace {
        medium,
        area_m2,
        pairs,
        gas_context: context,
        solve: DustyGasFaceSolve {
            maximum_temperature_difference_k: 1e-6,
            maximum_relative_pressure_difference: 1e-6,
            gradient_tolerance_per_m: 1e-7,
        },
    }
}
fn energy() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-8,
        maximum_iterations: 128,
    }
}
fn controls(dt: f64) -> PoreDomainSolve {
    PoreDomainSolve {
        maximum_step_s: dt,
        minimum_step_s: 1e-9,
        maximum_substeps: 512,
        maximum_attempts: 2048,
        maximum_throat_volume_fraction: 0.01,
        transport: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.1,
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 1e-7,
            },
            entropy_tolerance_j_k: 1e-10,
        },
    }
}
fn fixture() -> (PoreDomain, [DustyGasMedium; 2], MaterialCard, QueryPoint) {
    let cell = PoreCellGeometry::try_new(medium(1e-6, 0.4, 2.0, 0.025)).unwrap();
    let volume = cell.volume_m3(BODY_AREA).unwrap();
    let pore = cell
        .initialize_gas(BODY_AREA, gas_model(0.2), point(T, P))
        .unwrap();
    let domain = PoreDomain::try_new(
        BODY_AREA,
        vec![cell],
        gas(0.8, volume),
        vec![pore],
        gas(0.2, volume),
    )
    .unwrap();
    (
        domain,
        [medium(1e-6, 0.4, 2.0, 0.001), medium(1e-6, 0.4, 2.0, 0.001)],
        coefficient(),
        point(T, P),
    )
}
fn equivalent_diffusivity(medium: &DustyGasMedium) -> f64 {
    // Independent binary equal-mass DGM reduction: resistances add.
    let factor = medium.effective_diffusivity_factor();
    let dk = factor
        * (2.0 * medium.mean_gas_pore_radius_m() / 3.0)
        * (8.0 * R * T / (core::f64::consts::PI * MASS)).sqrt();
    1.0 / (1.0 / dk + 1.0 / (factor * D))
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn pore_inventory_delays_forwarding_and_carries_species_enthalpy_once() {
    let (mut domain, throats, card, context) = fixture();
    let original = domain.clone();
    let pairs = [binding(&card, &context)];
    let faces = throats
        .each_ref()
        .map(|m| face(m, THROAT_AREA, &pairs, &context));
    let dt = 0.1;
    let preview = domain
        .step_candidate(dt, &faces, &[energy(); 3], controls(dt), || false)
        .unwrap();
    assert_eq!(domain, original);
    assert!(x(preview.domain().pore_gas(0).unwrap()) > 0.2);
    assert_eq!(
        preview.domain().right_reservoir(),
        original.right_reservoir(),
        "arrivals cannot pass through a second throat in the original snapshot"
    );
    let q = (P / (R * T)) * MASS * equivalent_diffusivity(&throats[0]) * THROAT_AREA
        / throats[0].thickness_m()
        * (0.8 - 0.2);
    close(
        a_mass(preview.domain().pore_gas(0).unwrap()) - a_mass(original.pore_gas(0).unwrap()),
        q * dt,
        2e-14,
    );
    close(preview.step().neglected_throat_volume_m3, 8e-6, 1e-20);
    close(preview.step().throat_volume_fraction, 0.008, 1e-16);
    assert!(preview.step().entropy_change_j_k > 0.0);
    assert!(
        preview.step().total_energy_residual_j.abs()
            <= controls(dt).transport.exchange.energy_tolerance_j
    );
    let expected = preview.domain().clone();
    domain
        .advance(dt, &faces, &[energy(); 3], controls(dt), || false)
        .unwrap();
    assert_eq!(domain, expected);
    domain
        .advance(dt, &faces, &[energy(); 3], controls(dt), || false)
        .unwrap();
    assert!(x(domain.right_reservoir()) > 0.2);
    for id in 0..3 {
        close(
            domain.parcels()[&id].mass_kg(),
            original.parcels()[&id].mass_kg(),
            2e-14,
        );
        assert_eq!(
            domain.parcels()[&id].volume_m3(),
            original.parcels()[&id].volume_m3()
        );
    }
}

#[test]
fn fixed_step_refinement_converges_to_independent_three_volume_solution() {
    let (domain, throats, card, context) = fixture();
    let pairs = [binding(&card, &context)];
    let faces = throats
        .each_ref()
        .map(|m| face(m, THROAT_AREA, &pairs, &context));
    let duration = 4.0;
    let rate = equivalent_diffusivity(&throats[0]) * THROAT_AREA
        / (throats[0].thickness_m() * domain.pore_gas(0).unwrap().volume_m3());
    // Eigenvalues of the physical three-node path Laplacian are 0, rate, 3*rate.
    let (slow, fast) = ((-rate * duration).exp(), (-3.0 * rate * duration).exp());
    let exact = [
        0.4 + 0.3 * slow + 0.1 * fast,
        0.4 - 0.2 * fast,
        0.4 - 0.3 * slow + 0.1 * fast,
    ];
    let mut errors = Vec::new();
    for dt in [0.5, 0.25, 0.125] {
        let candidate = domain
            .step_candidate(duration, &faces, &[energy(); 3], controls(dt), || false)
            .unwrap();
        let error = (0..3)
            .map(|i| (x(&candidate.domain().parcels()[&(i as u64)]) - exact[i]).abs())
            .fold(0.0_f64, f64::max);
        errors.push(error);
        assert_eq!(candidate.domain().time_s(), duration);
        assert_eq!(
            candidate.step().attempts,
            candidate.step().accepted_substeps
        );
        assert!(
            candidate.step().maximum_species_mass_residual_kg
                <= controls(dt).transport.exchange.mass_tolerance_kg
        );
        assert!(candidate.step().entropy_change_j_k > 0.0);
    }
    assert!(
        errors[2] > 1e-7,
        "the oracle must resolve the first-order integration error"
    );
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.2,
            "{errors:?}"
        );
    }
}

#[test]
fn heterogeneous_body_storage_and_throats_follow_their_actual_source_coefficients() {
    let cells = [medium(2e-7, 0.3, 3.0, 0.02), medium(1e-6, 0.6, 2.0, 0.04)]
        .map(|m| PoreCellGeometry::try_new(m).unwrap());
    let pore = cells
        .iter()
        .zip([0.6, 0.3])
        .map(|(c, composition)| {
            c.initialize_gas(BODY_AREA, gas_model(composition), point(T, P))
                .unwrap()
        })
        .collect();
    let domain = PoreDomain::try_new(
        BODY_AREA,
        cells.to_vec(),
        gas(0.8, 0.001),
        pore,
        gas(0.2, 0.001),
    )
    .unwrap();
    let throats = [
        medium(1e-7, 0.3, 3.0, 0.001),
        medium(3e-7, 0.5, 4.0, 0.005),
        medium(2e-6, 0.6, 2.0, 0.002),
    ];
    let areas = [0.001, 0.002, 0.0005];
    let (card, context) = (coefficient(), point(T, P));
    let pairs = [binding(&card, &context)];
    let faces: Vec<_> = throats
        .iter()
        .zip(areas)
        .map(|(m, area)| face(m, area, &pairs, &context))
        .collect();
    let dt = 0.1;
    let result = domain
        .step_candidate(dt, &faces, &[energy(); 4], controls(dt), || false)
        .unwrap();
    let initial = [0.8, 0.6, 0.3, 0.2];
    let rates: Vec<_> = throats
        .iter()
        .enumerate()
        .map(|(i, m)| {
            P / (R * T) * MASS * equivalent_diffusivity(m) * areas[i] / m.thickness_m()
                * (initial[i] - initial[i + 1])
        })
        .collect();
    for i in 0..3 {
        let row = result
            .step()
            .last_substep
            .dusty_gas
            .iter()
            .find(|r| r.id == i as u64)
            .unwrap();
        assert_eq!(row.flux.medium, throats[i]);
        let flux = row
            .flux
            .species
            .iter()
            .find(|r| r.species == SpeciesId::new("A").unwrap())
            .unwrap();
        close(flux.mass_rate_left_to_right_kg_s, rates[i], 1e-16);
    }
    for i in 0..2 {
        let actual = result.domain().pore_gas(i).unwrap();
        close(
            actual.volume_m3(),
            cells[i].gas_filled_porosity() * BODY_AREA * cells[i].thickness_m(),
            1e-18,
        );
        close(
            a_mass(actual) - a_mass(domain.pore_gas(i).unwrap()),
            (rates[i] - rates[i + 1]) * dt,
            2e-14,
        );
    }
    let a = SpeciesId::new("A").unwrap();
    let mass_in = result
        .step()
        .faces
        .iter()
        .find(|r| r.face == 0 && r.species == a)
        .unwrap();
    let mass_out = result
        .step()
        .faces
        .iter()
        .find(|r| r.face == 2 && r.species == a)
        .unwrap();
    close(
        (mass_in.mass_left_to_right_kg - mass_in.mass_right_to_left_kg)
            - (mass_out.mass_left_to_right_kg - mass_out.mass_right_to_left_kg),
        (0..2)
            .map(|i| {
                a_mass(result.domain().pore_gas(i).unwrap()) - a_mass(domain.pore_gas(i).unwrap())
            })
            .sum(),
        2e-14,
    );
}

#[test]
fn adaptive_withdrawal_retries_preserve_window_atomicity_and_cancellation() {
    let (mut domain, throats, card, context) = fixture();
    let original = domain.clone();
    let pairs = [binding(&card, &context)];
    let faces = throats
        .each_ref()
        .map(|m| face(m, THROAT_AREA, &pairs, &context));
    let mut solve = controls(0.25);
    solve.transport.exchange.maximum_species_fraction = 0.002;
    let candidate = domain
        .step_candidate(0.25, &faces, &[energy(); 3], solve, || false)
        .unwrap();
    assert!(candidate.step().attempts > candidate.step().accepted_substeps);
    assert!(candidate.step().minimum_accepted_step_s < solve.maximum_step_s);
    let mut exhausted = solve;
    exhausted.maximum_substeps = 1;
    assert!(
        matches!(domain.advance(0.25, &faces, &[energy(); 3], exhausted, || false),
        Err(PoreDomainError::WorkLimit { elapsed_s, accepted_substeps: 1, .. }) if elapsed_s > 0.0)
    );
    assert_eq!(domain, original);
    let mut polls = 0;
    domain
        .clone()
        .advance(0.25, &faces, &[energy(); 3], solve, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert_eq!(
        domain.advance(0.25, &faces, &[energy(); 3], solve, || {
            count += 1;
            count == polls
        }),
        Err(PoreDomainError::Cancelled)
    );
    assert_eq!(domain, original);
}

#[test]
fn storage_geometry_and_neglected_throat_volume_are_explicit_refusals() {
    let (mut domain, throats, card, context) = fixture();
    let original = domain.clone();
    assert!(
        PoreDomain::try_new(
            BODY_AREA,
            domain.cells().to_vec(),
            domain.left_reservoir().clone(),
            vec![gas(0.2, 2.0 * domain.pore_gas(0).unwrap().volume_m3())],
            domain.right_reservoir().clone()
        )
        .is_err()
    );
    assert!(domain.cells()[0].volume_m3(0.0).is_err());
    let pairs = [binding(&card, &context)];
    let faces = throats
        .each_ref()
        .map(|m| face(m, THROAT_AREA, &pairs, &context));
    let mut solve = controls(0.1);
    solve.maximum_throat_volume_fraction = 0.001;
    assert!(matches!(
        domain.advance(0.1, &faces, &[energy(); 3], solve, || false),
        Err(PoreDomainError::Invalid(
            "neglected throat storage exceeds its declared pore-body volume fraction"
        ))
    ));
    assert_eq!(domain, original);
    // Closed throats have no active material queries, flux or omitted gas volume.
    let closed = throats.each_ref().map(|m| face(m, 0.0, &[], &context));
    solve.maximum_throat_volume_fraction = 0.0;
    let step = domain
        .advance(0.1, &closed, &[energy(); 3], solve, || false)
        .unwrap();
    assert_eq!(domain.parcels(), original.parcels());
    assert_eq!(step.neglected_throat_volume_m3, 0.0);
    assert_eq!(step.throat_volume_fraction, 0.0);
    assert!(step.faces.is_empty());
}
