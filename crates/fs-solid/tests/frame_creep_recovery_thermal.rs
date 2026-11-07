//! Prescribed thermal histories in actual recovering creep/yield members.
//! Independent scalar references are manufactured, not calibrated alloy data.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::recovery::plastic::source::RecoveringNortonJ2Source;
use fs_material::creep::recovery::source::{
    RecoveringNortonSource, recovering_norton_requirements,
};
use fs_material::creep::source::NortonSourceQuery;
use fs_material::creep::thermal::expansion::ExpandingCreepSource;
use fs_material::creep::{CreepTolerances, GAS_CONSTANT_J_MOL_K};
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    YIELD_STRESS_PROPERTY, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::fiber::material::FiberMaterial;
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::creep::FrameCreepTolerances;
use fs_solid::{Frame2d, Geometry, SectionModel};

const E: f64 = 1e9;
const STRESS: f64 = 1e6;
const RATE: f64 = 1e-5;
const PRIMARY: f64 = 4e-4;
const RECOVERY: f64 = 0.2;
const FLOW_Q: f64 = 20_000.0;
const RECOVERY_Q: f64 = 40_000.0;
const HARDENING: f64 = 4e7;
const EXPANSION: f64 = 2e-5;
const AREA: f64 = 5e-4;
const VOLUME: f64 = 2.0 * AREA;
const ARM: f64 = 0.05;

fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_spec(), t)
        .unwrap()
}

fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            temperature_spec(),
            300.0,
            500.0,
        ),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured thermal recovering frame".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (requirement, value) in recovering_norton_requirements().unwrap().iter().zip([
        E, 0.25, STRESS, RATE, 1.0, FLOW_Q, 300.0, PRIMARY, 1.0, 0.0, RECOVERY, RECOVERY_Q, 300.0,
    ]) {
        claims
            .insert_claim(claim(requirement.name(), requirement.dims(), value))
            .unwrap();
    }
    let mut yield_claim = claim(YIELD_STRESS_PROPERTY, Pressure::DIMS, 8e5);
    yield_claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 8e5), (500.0, 4e5)],
        dims: Pressure::DIMS,
    };
    yield_claim.interpolation = InterpolationPolicy::LinearInside;
    for property in [
        yield_claim,
        claim(
            ISOTROPIC_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            HARDENING,
        ),
        claim(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            EXPANSION,
        ),
    ] {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal recovering member".into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn sources(
    card: &MaterialCard,
) -> (
    ExpandingCreepSource<RecoveringNortonSource<'_>>,
    ExpandingCreepSource<RecoveringNortonJ2Source<'_>>,
) {
    let lower = point(300.0);
    let upper = point(500.0);
    let query = || NortonSourceQuery {
        lower: &lower,
        upper: &upper,
        temperature_axis: "T",
        selection: MaterialPropertySelection::SingleClaimOnly,
        strain_limit: 0.03,
    };
    let pure = RecoveringNortonSource::new(card, query(), || false).unwrap();
    let joint = RecoveringNortonJ2Source::new(
        card,
        query(),
        HardeningPropertyBasis::PlasticStrainModulus,
        || false,
    )
    .unwrap();
    let expansion = resolve_material_state_point(
        card,
        &lower,
        &[ScalarPropertyRequirement::try_new(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (
        ExpandingCreepSource::new(card, pure, &expansion, || false).unwrap(),
        ExpandingCreepSource::new(card, joint, &expansion, || false).unwrap(),
    )
}

fn frame(law: FiberLaw, axial_only: bool) -> Frame2d {
    let material = FiberMaterial::authored(law, 1000.0).unwrap();
    assert!(material.source().is_none());
    let mut frame = Frame2d::new(Geometry::Linear);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    if axial_only {
        frame.fix(1, [false, true, true]);
    }
    frame
        .add_element(0, 1, 3, 1000.0 * VOLUME, &|| {
            SectionModel::Fiber(Section::new(
                [-ARM, ARM]
                    .into_iter()
                    .map(|y| Fiber {
                        y,
                        area: AREA,
                        law: material.law().clone(),
                    })
                    .collect(),
            ))
        })
        .unwrap();
    frame
}

fn controls() -> FrameCreepTolerances {
    FrameCreepTolerances {
        material: CreepTolerances {
            stress_pa: 1e-3,
            energy_j_m3: 1e-5,
        },
        energy_j: 1e-7,
        force_n: 1e-5,
        moment_n_m: 1e-6,
        maximum_fiber_points: 6,
    }
}

fn load(stress: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, VOLUME * stress, 0.0, 0.0]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

fn recovery_rate(t: f64) -> f64 {
    RECOVERY * (RECOVERY_Q / GAS_CONSTANT_J_MOL_K * (1.0 / 300.0 - 1.0 / t)).exp()
}

fn yield_stress(t: f64) -> f64 {
    8e5 - 2e3 * (t - 300.0)
}

fn flow(stress: f64, resistance: f64, t: f64, dt: f64) -> f64 {
    let a = PRIMARY * (1.0 + recovery_rate(t) * dt);
    let rate = RATE * stress.abs() / STRESS
        * (FLOW_Q / GAS_CONSTANT_J_MOL_K * (1.0 / 300.0 - 1.0 / t)).exp();
    let b = a + resistance;
    let c = dt * rate * a;
    2.0 * c / (b + (b * b + 4.0 * c).sqrt())
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct History {
    creep: f64,
    alpha_c: f64,
    plastic: f64,
    alpha_p: f64,
    resistance: f64,
    temperature: f64,
}

fn observe(
    frame: &Frame2d,
    joint: bool,
    section: usize,
    fiber: usize,
) -> (History, &ResolvedMaterialStatePoint) {
    if joint {
        let (law, state) = frame
            .thermal_creep_recovery_plastic_fiber(0, section, fiber)
            .unwrap();
        let m = &state.mechanical;
        (
            History {
                creep: m.creep_strain(),
                alpha_c: m.equivalent_creep_strain(),
                plastic: m.plastic_strain(),
                alpha_p: m.equivalent_plastic_strain(),
                resistance: m.resistance_strain(),
                temperature: state.temperature_k(),
            },
            law.material().unwrap(),
        )
    } else {
        let (law, state) = frame
            .thermal_creep_recovery_fiber(0, section, fiber)
            .unwrap();
        let m = &state.mechanical;
        (
            History {
                creep: m.creep_strain(),
                alpha_c: m.equivalent_creep_strain(),
                plastic: 0.0,
                alpha_p: 0.0,
                resistance: m.resistance_strain(),
                temperature: state.temperature_k(),
            },
            law.material().unwrap(),
        )
    }
}

#[test]
fn thermal_members_preserve_three_histories_and_reciprocal_heat_through_dwell_reload() {
    let card = card();
    let (pure, combined) = sources(&card);
    for joint in [false, true] {
        let at = |t, dt| {
            if joint {
                FiberLaw::ThermalRecoveringNortonJ2(
                    combined.at_uniaxial_temperature(t, dt).unwrap(),
                )
            } else {
                FiberLaw::ThermalRecoveringNorton(pure.at_uniaxial_temperature(t, dt).unwrap())
            }
        };
        let mut accepted = frame(at(300.0, 0.5), true);
        let receipts: Vec<_> = observe(&accepted, joint, 0, 0)
            .1
            .properties()
            .iter()
            .map(|p| p.answer().receipt.selected)
            .collect();
        assert_eq!(receipts.len(), if joint { 15 } else { 13 });
        let (mut creep, mut plastic, mut resistance) = (0.0, 0.0_f64, 0.0);
        let (mut old_t, mut old_stress) = (300.0, 0.0);
        for (index, (t, dt, stress)) in [
            (300.0, 0.5, 1e6),
            (350.0, 0.5, 1e6),
            (350.0, 1.0, 0.0),
            (400.0, 2.0, 0.0),
            (400.0, 0.5, 1.1e6),
        ]
        .into_iter()
        .enumerate()
        {
            let before = observe(&accepted, joint, 0, 0).0;
            let prepared = accepted
                .with_creep_laws(&vec![at(t, dt); 6], "T", || false)
                .unwrap();
            assert_eq!(observe(&prepared, joint, 0, 0).0, before);
            assert_eq!(prepared.displacements(), accepted.displacements());
            let dc = flow(stress, resistance, t, dt);
            let dp = if joint {
                ((stress - yield_stress(t)) / HARDENING - plastic).max(0.0)
            } else {
                0.0
            };
            let next_r = (resistance + dc) / (1.0 + recovery_rate(t) * dt);
            let step = prepared
                .solve_creep_step(&load(stress), controls(), || false)
                .unwrap();
            assert_eq!(observe(&prepared, joint, 0, 0).0, before);
            assert_eq!(observe(&accepted, joint, 0, 0).0, before);
            assert_eq!(step.duration_s, dt);
            close(
                step.frame.displacements()[3],
                stress / E + creep + dc + plastic + dp + EXPANSION * (t - 300.0),
                3e-11,
            );
            close(step.frame.element_forces(0)[0], VOLUME * stress, 2e-5);
            for section in 0..3 {
                for fiber in 0..2 {
                    let (h, material) = observe(&step.frame, joint, section, fiber);
                    close(h.creep, creep + dc, 3e-11);
                    close(h.alpha_c, creep + dc, 3e-11);
                    close(h.plastic, plastic + dp, 3e-11);
                    close(h.alpha_p, plastic + dp, 3e-11);
                    close(h.resistance, next_r, 3e-12);
                    assert_eq!(h.temperature, t);
                    assert_eq!(material.state_coordinate("T"), Some(t));
                    assert_eq!(material.card_identity(), card.content_hash());
                    assert_eq!(
                        material
                            .properties()
                            .iter()
                            .map(|p| p.answer().receipt.selected)
                            .collect::<Vec<_>>(),
                        receipts
                    );
                }
            }
            let energy = step.energy;
            close(energy.creep_dissipation_j, VOLUME * stress * dc, 2e-8);
            close(
                energy.plastic_dissipation_j,
                VOLUME * yield_stress(t) * dp,
                2e-8,
            );
            close(
                energy.physical_dissipation_j,
                VOLUME * (stress * dc + yield_stress(t) * dp),
                2e-8,
            );
            close(
                energy.final_stored_j,
                VOLUME * (0.5 * stress * stress / E + 0.5 * HARDENING * (plastic + dp).powi(2)),
                2e-8,
            );
            close(
                energy.reversible_heat_j,
                -VOLUME * old_t * EXPANSION * (stress - old_stress),
                2e-8,
            );
            close(
                energy.thermoelastic_storage_change_j,
                VOLUME * EXPANSION * (t * stress - old_t * old_stress),
                2e-8,
            );
            close(
                energy.caloric_heat_j(),
                energy.physical_dissipation_j + energy.reversible_heat_j,
                0.0,
            );
            close(energy.balance_residual_j, 0.0, controls().energy_j);
            if stress == 0.0 {
                assert!(next_r < resistance);
                close(energy.physical_dissipation_j, 0.0, 2e-10);
                if index == 3 {
                    close(energy.caloric_heat_j(), 0.0, 2e-8);
                    close(energy.applied_work_j + energy.support_work_j, 0.0, 2e-8);
                }
            } else if joint {
                assert!(dp > 0.0 && dc > 0.0);
            }
            (creep, plastic, resistance, old_t, old_stress) =
                (creep + dc, plastic + dp, next_r, t, stress);
            accepted = step.frame;
        }
    }
}

#[test]
fn unloaded_temperature_gradient_bends_a_crept_member_and_recovers_each_fiber_locally() {
    let card = card();
    let (source, _) = sources(&card);
    let cold = source.at_uniaxial_temperature(300.0, 0.5).unwrap();
    let loaded = frame(FiberLaw::ThermalRecoveringNorton(cold), false)
        .solve_creep_step(&load(STRESS), controls(), || false)
        .unwrap()
        .frame;
    let released = loaded
        .solve_creep_step(&load(0.0), controls(), || false)
        .unwrap()
        .frame;
    let before = observe(&released, false, 0, 0).0;
    assert!(before.creep > 0.0 && before.resistance > 0.0);
    let bottom = source.at_uniaxial_temperature(300.0, 2.0).unwrap();
    let top = source.at_uniaxial_temperature(400.0, 2.0).unwrap();
    let laws: Vec<_> = (0..3).flat_map(|_| [bottom.clone(), top.clone()]).collect();
    let prepared = released.with_creep_laws(&laws, "T", || false).unwrap();
    let step = prepared
        .solve_creep_step(&[0.0; 6], controls(), || false)
        .unwrap();
    let free = EXPANSION * 100.0;
    let curvature = -free / (2.0 * ARM);
    close(
        step.frame.displacements()[3],
        before.creep + 0.5 * free,
        4e-11,
    );
    close(step.frame.displacements()[4], 0.5 * curvature, 4e-10);
    close(step.frame.displacements()[5], curvature, 4e-10);
    for section in 0..3 {
        for (fiber, t) in [300.0, 400.0].into_iter().enumerate() {
            let (h, _) = observe(&step.frame, false, section, fiber);
            close(h.creep, before.creep, 3e-11);
            close(h.alpha_c, before.alpha_c, 3e-11);
            close(
                h.resistance,
                before.resistance / (1.0 + 2.0 * recovery_rate(t)),
                3e-12,
            );
            assert_eq!(h.temperature, t);
        }
    }
    assert!(
        observe(&step.frame, false, 0, 1).0.resistance
            < observe(&step.frame, false, 0, 0).0.resistance
    );
    for force in step.frame.element_forces(0) {
        close(force, 0.0, 2e-5);
    }
    close(step.energy.physical_dissipation_j, 0.0, 2e-10);
    close(step.energy.reversible_heat_j, 0.0, 2e-8);
    close(step.energy.balance_residual_j, 0.0, controls().energy_j);
    assert_eq!(observe(&released, false, 0, 0).0, before);
}
