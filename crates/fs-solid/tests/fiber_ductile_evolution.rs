//! G0/G1/G3 native section response and original source handoff.
//! Synthetic coefficients test implementation, not a measured alloy.

use fs_material::IsotropicElastic;
use fs_material::ductile::evolution::{
    DuctileJ2Conditions, DuctileJ2Settings, DuctileWorkEvolution,
};
use fs_material::ductile::uniaxial::{UniaxialDuctileJ2, UniaxialDuctileJ2State};
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::UniaxialJ2;
use fs_solid::fiber::ductile::DuctileFiberLaw;
use fs_solid::fiber::material::FiberMaterial;
use fs_solid::fiber::{Fiber, FiberLaw, FiberState, Section};

const E: f64 = 1000.0;
const ARM: f64 = 0.2;
const AREA: f64 = 0.5;

fn native() -> UniaxialJ2 {
    UniaxialJ2::new(
        J2Plasticity::new(IsotropicElastic::new(E, 0.3, 0.05).unwrap(), 1.0, 100.0).unwrap(),
    )
    .unwrap()
}
fn settings() -> DuctileJ2Settings {
    DuctileJ2Settings {
        characteristic_length_m: 0.1,
        post_initiation_plastic_work_j_m2: 0.001,
        maximum_damage_increment: 1.0,
        // EnergyTolerances has no Default (the ductile settings take no material
        // defaults); these are frame_ductile_damage.rs's, the sibling test of the
        // same commit on the same UniaxialDuctileJ2 -> DuctileFiberLaw path.
        energy_tolerances: EnergyTolerances {
            energy_j_m3: 1e-9,
            yield_pa: 1e-8,
        },
    }
}
fn conditions() -> DuctileJ2Conditions {
    DuctileJ2Conditions {
        duration_s: 0.1,
        temperature_k: 300.0,
    }
}
fn law(profile: DuctileWorkEvolution, controls: DuctileJ2Settings) -> DuctileFiberLaw {
    let initiation = TabulatedDuctileLocus::new(
        vec![(-0.5, 0.002), (0.5, 0.002)],
        DuctileDomain {
            triaxiality: [-0.5, 0.5],
            strain_rate_s_inv: [1e-12, 1e8],
            temperature_k: [300.0, 300.0],
        },
    )
    .unwrap();
    UniaxialDuctileJ2::new(native(), initiation, controls)
        .unwrap()
        .with_evolution(profile)
        .prepared(conditions())
        .unwrap()
        .into()
}
fn make_section(law: &DuctileFiberLaw) -> Section {
    Section::new(
        [-ARM, ARM]
            .map(|y| Fiber {
                y,
                area: AREA,
                law: FiberLaw::DuctileJ2(law.clone()),
            })
            .to_vec(),
    )
}
fn state(section: &Section, index: usize) -> &UniaxialDuctileJ2State {
    let FiberState::DuctileJ2(state) = &section.states[index] else {
        panic!("native history")
    };
    state
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn native_damage_enters_section_force_and_consistent_tangent_without_flooring() {
    for profile in [
        DuctileWorkEvolution::Exponential,
        DuctileWorkEvolution::FiniteSupportQuadratic,
    ] {
        let law = law(profile, settings());
        let mut section = make_section(&law);
        let strain = 0.008;
        let trial = law.response(strain, state(&section, 0)).unwrap();
        assert!(trial.state.damage() > 0.0 && trial.tangent_pa < 0.0);
        let response = section.try_respond(strain, 0.0).unwrap();
        close(response.n, 2.0 * AREA * trial.stress_pa, 1e-13);
        close(response.m, 0.0, 1e-13);
        close(response.tangent[0][0], 2.0 * AREA * trial.tangent_pa, 1e-13);
        close(
            response.tangent[1][1],
            2.0 * AREA * ARM * ARM * trial.tangent_pa,
            1e-13,
        );
        close(response.tangent[0][1], 0.0, 1e-13);
        // Differentiate the actual assembled section response at frozen old history.
        let delta = 1e-8;
        for column in 0..2 {
            let (ep, kp, em, km) = if column == 0 {
                (strain + delta, 0.0, strain - delta, 0.0)
            } else {
                (strain, delta, strain, -delta)
            };
            let plus = section.try_respond(ep, kp).unwrap();
            let minus = section.try_respond(em, km).unwrap();
            close(
                (plus.n - minus.n) / (2.0 * delta),
                response.tangent[0][column],
                2e-7,
            );
            close(
                (plus.m - minus.m) / (2.0 * delta),
                response.tangent[1][column],
                2e-7,
            );
        }
        assert_eq!(
            state(&section, 0).damage(),
            0.0,
            "a response must not commit history"
        );
        section.try_commit(strain, 0.0).unwrap();
        assert_eq!(state(&section, 0), &trial.state);
        assert_eq!(state(&section, 1), &trial.state);
        let damaged = state(&section, 0).clone();
        let unload = section.try_respond(strain - 0.0002, 0.0).unwrap();
        close(
            unload.tangent[0][0],
            2.0 * AREA * E * damaged.integrity(),
            1e-12,
        );
        section.try_commit(strain - 0.0002, 0.0).unwrap();
        assert_eq!(state(&section, 0).damage(), damaged.damage());
        assert_eq!(
            state(&section, 0).equivalent_plastic_strain(),
            damaged.equivalent_plastic_strain()
        );
        section.try_commit(-strain, 0.0).unwrap();
        assert!(state(&section, 0).damage() > damaged.damage());
        assert!(
            state(&section, 0).equivalent_plastic_strain() > damaged.equivalent_plastic_strain()
        );
        assert!(
            state(&section, 0).last_step_energy().residual_j_m3.abs()
                <= settings().energy_tolerances.energy_j_m3
        );
    }
}

#[test]
fn finite_support_completion_retains_history_and_adds_no_residual_section_stiffness() {
    let controls = DuctileJ2Settings {
        post_initiation_plastic_work_j_m2: 0.0001,
        ..settings()
    };
    let law = law(DuctileWorkEvolution::FiniteSupportQuadratic, controls);
    let mut section = make_section(&law);
    section.try_commit(0.008, 0.0).unwrap();
    let failed = state(&section, 0).clone();
    assert!(failed.is_fully_failed());
    let alpha = failed.equivalent_plastic_strain();
    for strain in [0.008, 0.007, -0.008] {
        let response = section.try_respond(strain, 0.0).unwrap();
        assert_eq!(response.n, 0.0);
        assert_eq!(response.m, 0.0);
        assert_eq!(response.tangent, [[0.0; 2]; 2]);
        section.try_commit(strain, 0.0).unwrap();
        assert!(state(&section, 0).is_fully_failed());
        assert!(state(&section, 0).equivalent_plastic_strain() >= alpha);
        assert_eq!(section.fibers.len(), 2);
        assert_eq!(section.fibers[0].area, AREA);
    }
    // Force-based member equilibrium must diagnose the resulting mechanism;
    // the material owner neither deletes fibers nor manufactures stiffness.
    assert_eq!(
        state(&section, 0).degraded_equivalent_plastic_work_j_m2(),
        controls.post_initiation_plastic_work_j_m2
    );
}

#[test]
fn interval_rebinding_preserves_identity_and_late_refusal_preserves_the_whole_section() {
    let first = law(DuctileWorkEvolution::Exponential, settings());
    let faster = first.with_duration(0.01).unwrap();
    assert!(first.same_model(&faster));
    assert!(faster.owns_state(&first.initial_state()));
    assert_eq!(
        faster.conditions().temperature_k,
        first.conditions().temperature_k
    );
    let independent = law(DuctileWorkEvolution::Exponential, settings());
    assert!(!first.same_model(&independent));
    let mut section = make_section(&first);
    section.states[1] = FiberState::DuctileJ2(independent.initial_state());
    let old = [state(&section, 0).clone(), state(&section, 1).clone()];
    assert!(section.try_commit(0.008, 0.0).is_err());
    assert_eq!(state(&section, 0), &old[0]);
    assert_eq!(state(&section, 1), &old[1]);
    let restricted = law(
        DuctileWorkEvolution::Exponential,
        DuctileJ2Settings {
            maximum_damage_increment: 0.01,
            ..settings()
        },
    );
    let mut section = make_section(&restricted);
    assert!(section.try_commit(0.008, 0.0).is_err());
    assert_eq!(state(&section, 0).damage(), 0.0);
    assert_eq!(state(&section, 1).damage(), 0.0);
}

mod source {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId,
        LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::ductile::source::{
        JOHNSON_COOK_INITIATION_LAW_ID, PLASTIC_STRAIN_RATE_AXIS, STRESS_TRIAXIALITY_AXIS,
        resolve_johnson_cook_initiation,
    };
    use fs_material::plastic::source::HardeningPropertyBasis;
    use fs_material::state_point::MaterialPropertySelection;
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};

    fn temperature() -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ))
    }
    fn point(value: f64) -> QueryPoint {
        QueryPoint::new()
            .with_quantity("T", temperature(), value)
            .unwrap()
    }
    fn material_card(revision: u32) -> MaterialCard {
        let fixed = ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 600.0);
        let artifact = fs_blake3::hash_bytes(b"synthetic progressive fiber JC fit");
        let provenance = Provenance {
            source: "synthetic native progressive fiber".into(),
            license: "test-only".into(),
            artifact: Some(artifact),
        };
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("density", Density::DIMS, 7800.0),
            ("young_modulus", Pressure::DIMS, E),
            ("poisson_ratio", Dims::NONE, 0.3),
            ("yield_stress", Pressure::DIMS, 1.0),
            ("isotropic_hardening_modulus", Pressure::DIMS, 100.0),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: fixed.clone(),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: provenance.clone(),
                })
                .unwrap();
        }
        let rate_dims = Dims([0, 0, -1, 0, 0, 0]);
        let model = ConstitutiveModelCard {
            law: LawId(JOHNSON_COOK_INITIATION_LAW_ID.into()),
            law_version: 1,
            parameters: [
                ("d1", 0.002, Dims::NONE),
                ("d2", 0.0, Dims::NONE),
                ("d3", 0.0, Dims::NONE),
                ("d4", 0.0, Dims::NONE),
                ("d5", 0.0, Dims::NONE),
                ("reference_strain_rate_s_inv", 1.0, rate_dims),
                ("transition_temperature_k", 300.0, Temperature::DIMS),
                ("melting_temperature_k", 1000.0, Temperature::DIMS),
            ]
            .into_iter()
            .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
            state_schema_version: 1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: fixed
                .with_quantity(
                    STRESS_TRIAXIALITY_AXIS,
                    QuantitySpec::dimensional(Dims::NONE),
                    -0.5,
                    0.5,
                )
                .with_quantity(
                    PLASTIC_STRAIN_RATE_AXIS,
                    QuantitySpec::dimensional(rate_dims),
                    1e-12,
                    1e8,
                ),
            sources: vec![artifact],
            provenance,
        };
        // A lineage starts at revision 0; each later revision supersedes the one before.
        let mut card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic progressive frame alloy".into(),
                phase: "solid".into(),
                process: "synthetic fixed condition".into(),
                revision: 0,
            },
            claims.clone(),
            vec![model.clone()],
        )
        .unwrap();
        for _ in 0..revision {
            card = MaterialCard::supersede(&card, claims.clone(), vec![model.clone()]).unwrap();
        }
        card
    }
    fn resolve(card: &MaterialCard, point: &QueryPoint) -> FiberMaterial {
        FiberMaterial::resolve_j2(
            card,
            point,
            "density",
            HardeningPropertyBasis::PlasticStrainModulus,
            0.05,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    }

    #[test]
    fn original_mechanical_source_and_actual_temperature_survive_damage_admission() {
        let card = material_card(0);
        let cold = resolve(&card, &point(300.0));
        let hot = resolve(&card, &point(600.0));
        let query = point(300.0)
            .with_quantity(
                STRESS_TRIAXIALITY_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                1.0 / 3.0,
            )
            .unwrap()
            .with_quantity(
                PLASTIC_STRAIN_RATE_AXIS,
                QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                1.0,
            )
            .unwrap();
        let initiation =
            resolve_johnson_cook_initiation(&card, card.models()[0].content_hash(), &query, "T")
                .unwrap();
        let FiberLaw::J2(backbone) = cold.law() else {
            panic!("native J2 source")
        };
        let model = UniaxialDuctileJ2::new(backbone.clone(), initiation, settings()).unwrap();
        let prepared = model.prepared(conditions()).unwrap();
        let bound = cold
            .clone()
            .with_ductile_evolution(prepared.clone())
            .unwrap();
        assert_eq!(bound.density_kg_m3(), cold.density_kg_m3());
        assert_eq!(bound.source(), cold.source());
        let FiberLaw::DuctileJ2(native) = bound.law() else {
            panic!("progressive source")
        };
        assert!(
            native
                .response(0.008, &native.initial_state())
                .unwrap()
                .state
                .damage()
                > 0.0
        );
        let hot_law = model
            .prepared(DuctileJ2Conditions {
                temperature_k: 600.0,
                ..conditions()
            })
            .unwrap();
        assert!(
            cold.clone()
                .with_ductile_evolution(hot_law.clone())
                .is_err()
        );
        assert!(hot.with_ductile_evolution(hot_law).is_ok());
        let foreign = resolve(&material_card(1), &point(300.0));
        assert!(foreign.with_ductile_evolution(prepared.clone()).is_err());
        let pressure =
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
        let other_environment = resolve(
            &card,
            &point(300.0)
                .with_quantity("ambient-pressure", pressure, 101325.0)
                .unwrap(),
        );
        assert!(
            other_environment
                .with_ductile_evolution(prepared.clone())
                .is_err()
        );
        assert!(FiberMaterial::authored(FiberLaw::DuctileJ2(prepared.into()), 7800.0).is_err());
        assert!(
            cold.with_ductile_evolution(law(DuctileWorkEvolution::Exponential, settings()))
                .is_err()
        );
        // Gp and profile were supplied explicitly; no source property was fabricated.
        assert!(bound.source().unwrap().property(
            fs_material::ductile::evolution::source::POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY
        ).is_none());
    }
}
