//! Independent force/curvature references for native joint creep/plastic frames.
//! Manufactured coefficients establish implementation behavior, not service life.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::IsotropicElastic;
use fs_material::creep::plastic::NortonJ2;
use fs_material::creep::plastic::source::NortonJ2Source;
use fs_material::creep::plastic::uniaxial::UniaxialNortonJ2;
use fs_material::creep::source::{NortonSourceQuery, norton_requirements};
use fs_material::creep::{CreepTolerances, NortonCreep, NortonCreepParameters};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::state_point::{MaterialPropertySelection, YIELD_STRESS_PROPERTY};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::fiber::material::FiberMaterial;
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::creep::{FrameCreepError, FrameCreepTolerances};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel, SolidError};
use std::cell::Cell;

const E: f64 = 1e9;
const STRESS: f64 = 1e6;
const YIELD: f64 = 0.8e6;
const HARDENING: f64 = 2e8;
const RATE: f64 = 1e-5;
const DT: f64 = 10.0;
const AREA: f64 = 5e-4;
const ARM: f64 = 0.05;
const FORCE: f64 = 2.0 * AREA * STRESS;

fn law(dt: f64, temperature: f64) -> UniaxialNortonJ2 {
    let elastic = IsotropicElastic::new(E, 0.25, 0.05).unwrap();
    let creep = NortonCreep::new(
        elastic.clone(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: 1.0,
            activation_energy_j_mol: 20_000.0,
            reference_temperature_k: 300.0,
            temperature_range_k: [280.0, 360.0],
        },
        temperature,
        dt,
    )
    .unwrap();
    UniaxialNortonJ2::new(
        NortonJ2::new(creep, J2Plasticity::new(elastic, YIELD, HARDENING).unwrap()).unwrap(),
    )
}
fn frame(law: &UniaxialNortonJ2, geometry: Geometry, axial: bool) -> Frame2d {
    // Exercise the same explicit-density material dispatch used by native sections.
    let material = FiberMaterial::authored(FiberLaw::NortonJ2(law.clone()), 1000.0).unwrap();
    let mut frame = Frame2d::new(geometry);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    if axial {
        frame.fix(1, [false, true, true]);
    }
    frame
        .add_element(0, 1, 3, 1.0, &|| {
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
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}
fn axial_load(force: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, force, 0.0, 0.0]
}

#[test]
fn loaded_member_yields_creeps_relaxes_when_held_and_keeps_both_permanent_strains() {
    let mut accepted = frame(&law(DT, 300.0), Geometry::Linear, true);
    let plastic = (STRESS - YIELD) / HARDENING;
    for i in 1..=3 {
        let original = accepted.displacements().to_vec();
        let step = accepted
            .solve_creep_step(&axial_load(FORCE), controls(), || false)
            .unwrap();
        assert_eq!(accepted.displacements(), original);
        close(
            step.frame.displacements()[3],
            STRESS / E + plastic + i as f64 * DT * RATE,
            2e-11,
        );
        close(step.frame.element_forces(0)[0], FORCE, 1e-6);
        close(
            step.energy.creep_dissipation_j,
            2.0 * AREA * STRESS * DT * RATE,
            1e-8,
        );
        close(
            step.energy.plastic_dissipation_j,
            if i == 1 {
                2.0 * AREA * YIELD * plastic
            } else {
                0.0
            },
            1e-8,
        );
        close(
            step.energy.physical_dissipation_j,
            step.energy.creep_dissipation_j + step.energy.plastic_dissipation_j,
            1e-12,
        );
        close(
            step.energy.final_stored_j,
            2.0 * AREA * (0.5 * STRESS * STRESS / E + 0.5 * HARDENING * plastic * plastic),
            1e-8,
        );
        close(step.energy.balance_residual_j, 0.0, controls().energy_j);
        for p in 0..3 {
            for f in 0..2 {
                let (_, state) = step.frame.creep_plastic_fiber(0, p, f).unwrap();
                close(state.creep_strain(), i as f64 * DT * RATE, 2e-11);
                close(state.plastic_strain(), plastic, 2e-11);
                close(state.equivalent_creep_strain(), i as f64 * DT * RATE, 2e-11);
                close(state.equivalent_plastic_strain(), plastic, 2e-11);
            }
        }
        accepted = step.frame;
    }
    let mut held = accepted.clone();
    held.fix(1, [true; 3]);
    let fixed = held.displacements().to_vec();
    let ratio = 1.0 / (1.0 + E * DT * RATE / STRESS);
    for i in 1..=3 {
        let step = held
            .solve_creep_step(&axial_load(0.0), controls(), || false)
            .unwrap();
        assert_eq!(step.frame.displacements(), fixed);
        close(step.frame.element_forces(0)[0], FORCE * ratio.powi(i), 1e-6);
        close(
            step.frame
                .creep_plastic_fiber(0, 1, 0)
                .unwrap()
                .1
                .plastic_strain(),
            plastic,
            2e-11,
        );
        assert!(step.energy.creep_dissipation_j > 0.0);
        close(step.energy.plastic_dissipation_j, 0.0, 1e-8);
        close(
            step.energy.applied_work_j + step.energy.support_work_j,
            0.0,
            1e-14,
        );
        close(
            step.energy.initial_stored_j - step.energy.final_stored_j,
            step.energy.physical_dissipation_j + step.energy.numerical_dissipation_j,
            1e-8,
        );
        held = step.frame;
    }
    let released = accepted
        .solve_creep_step(&axial_load(0.0), controls(), || false)
        .unwrap();
    close(
        released.frame.displacements()[3],
        plastic + 3.0 * DT * RATE,
        2e-11,
    );
    close(released.frame.element_forces(0)[0], 0.0, 1e-6);
    close(released.energy.physical_dissipation_j, 0.0, 1e-8);
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn claim(name: &str, dims: Dims, lower: f64, upper: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, lower), (360.0, upper)],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            280.0,
            360.0,
        ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured joint frame fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn card(name: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (requirement, value) in norton_requirements()
        .unwrap()
        .iter()
        .zip([E, 0.25, STRESS, RATE, 1.0, 20_000.0, 300.0])
    {
        claims
            .insert_claim(claim(requirement.name(), requirement.dims(), value, value))
            .unwrap();
    }
    claims
        .insert_claim(claim(
            ISOTROPIC_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            HARDENING,
            HARDENING,
        ))
        .unwrap();
    claims
        .insert_claim(claim(YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.3e6, 0.5e6))
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn source(card: &MaterialCard) -> NortonJ2Source<'_> {
    NortonJ2Source::new(
        card,
        NortonSourceQuery {
            lower: &point(280.0),
            upper: &point(360.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.05,
        },
        HardeningPropertyBasis::PlasticStrainModulus,
        || false,
    )
    .unwrap()
}

#[test]
fn nonuniform_sourced_heating_changes_yield_curvature_and_actual_cantilever_shape() {
    let original = card("frame source A");
    let family = source(&original);
    let initial = UniaxialNortonJ2::from_source(family.at_temperature(300.0, DT).unwrap());
    let moment = FORCE * ARM;
    let load = [0.0, 0.0, 0.0, 0.0, 0.0, moment];
    let cold = frame(&initial, Geometry::Linear, false)
        .solve_creep_step(&load, controls(), || false)
        .unwrap()
        .frame;
    let temperatures = [300.0, 315.0, 330.0];
    let laws: Vec<_> = temperatures
        .into_iter()
        .flat_map(|t| {
            let law = UniaxialNortonJ2::from_source(family.at_temperature(t, DT).unwrap());
            [law.clone(), law]
        })
        .collect();
    let prepared = cold.with_creep_laws(&laws, "T", || false).unwrap();
    for p in 0..3 {
        for f in 0..2 {
            assert_eq!(
                prepared
                    .creep_plastic_fiber(0, p, f)
                    .unwrap()
                    .1
                    .tensor_history(),
                cold.creep_plastic_fiber(0, p, f)
                    .unwrap()
                    .1
                    .tensor_history()
            );
        }
    }
    let step = prepared
        .solve_creep_step(&load, controls(), || false)
        .unwrap();
    let mut expected_rotation = 0.0;
    let mut expected_deflection = 0.0;
    let mut creep_heat = 0.0;
    let mut plastic_heat = 0.0;
    for (p, ((t, weight), x)) in temperatures
        .into_iter()
        .zip([1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0])
        .zip([0.0, 0.5, 1.0])
        .enumerate()
    {
        let dc = DT * RATE * (20_000.0 / 8.314_462_618_153_24 * (1.0 / 300.0 - 1.0 / t)).exp();
        let yield_stress = 1.1e6 - 1e4 * (t - 300.0);
        let plastic = ((STRESS - yield_stress) / HARDENING).max(0.0);
        let curvature = (STRESS / E + DT * RATE + dc + plastic) / ARM;
        close(step.frame.section_deformations(0)[p].0, 0.0, 2e-11);
        close(step.frame.section_deformations(0)[p].1, curvature, 3e-10);
        expected_rotation += weight * curvature;
        expected_deflection += weight * (1.0 - x) * curvature;
        creep_heat += 2.0 * AREA * weight * STRESS * dc;
        plastic_heat += 2.0 * AREA * weight * yield_stress * plastic;
        let (law, bottom) = step.frame.creep_plastic_fiber(0, p, 0).unwrap();
        let (_, top) = step.frame.creep_plastic_fiber(0, p, 1).unwrap();
        close(bottom.creep_strain(), DT * RATE + dc, 2e-11);
        close(top.creep_strain(), -bottom.creep_strain(), 2e-11);
        close(bottom.plastic_strain(), plastic, 2e-11);
        close(top.plastic_strain(), -plastic, 2e-11);
        let material = law.material().unwrap();
        assert_eq!(material.card_identity(), original.content_hash());
        close(material.state_coordinate("T").unwrap(), t, 1e-12);
        for property in material.properties() {
            assert_eq!(
                property.answer().receipt.selected,
                initial
                    .material()
                    .unwrap()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
    close(step.frame.displacements()[4], expected_deflection, 3e-10);
    close(step.frame.displacements()[5], expected_rotation, 3e-10);
    close(step.energy.creep_dissipation_j, creep_heat, 2e-8);
    close(step.energy.plastic_dissipation_j, plastic_heat, 2e-8);
    close(step.energy.balance_residual_j, 0.0, controls().energy_j);
    let foreign = card("frame source B");
    let replacement =
        UniaxialNortonJ2::from_source(source(&foreign).at_temperature(330.0, DT).unwrap());
    assert!(
        cold.with_creep_laws(&vec![replacement; 6], "T", || false)
            .is_err()
    );
}

#[test]
fn unsupported_routes_mismatched_clocks_and_failed_trials_preserve_both_memories() {
    let material = law(DT, 300.0);
    let mut accepted = frame(&material, Geometry::Linear, true);
    let before = accepted.displacements().to_vec();
    let history = accepted
        .creep_plastic_fiber(0, 0, 0)
        .unwrap()
        .1
        .tensor_history();
    assert!(matches!(
        accepted.static_load(&axial_load(FORCE), 2),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert!(matches!(
        accepted.displacement_control(&axial_load(FORCE), 3, &[0.001]),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert!(matches!(
        accepted.modal_analysis(),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert!(matches!(
        accepted.newmark(
            &[0.0, 0.0],
            DT / 2.0,
            (1.0, 0.0),
            Rayleigh { a0: 0.0, a1: 0.0 }
        ),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert!(
        frame(&material, Geometry::Corotational, true)
            .solve_creep_step(&axial_load(FORCE), controls(), || false)
            .is_err()
    );
    let mut unequal = vec![material; 6];
    unequal[3] = law(DT / 2.0, 300.0);
    let prepared = accepted.with_creep_laws(&unequal, "T", || false).unwrap();
    assert!(
        prepared
            .solve_creep_step(&axial_load(FORCE), controls(), || false)
            .is_err()
    );
    assert!(
        accepted
            .solve_creep_step(&axial_load(100.0 * FORCE), controls(), || false)
            .is_err()
    );
    let polls = Cell::new(0usize);
    assert!(matches!(
        accepted.solve_creep_step(&axial_load(FORCE), controls(), || {
            polls.set(polls.get() + 1);
            polls.get() > 4
        }),
        Err(FrameCreepError::Cancelled)
    ));
    assert_eq!(accepted.displacements(), before);
    assert_eq!(
        accepted
            .creep_plastic_fiber(0, 0, 0)
            .unwrap()
            .1
            .tensor_history(),
        history
    );
    let next = accepted
        .solve_creep_step(&axial_load(FORCE), controls(), || false)
        .unwrap();
    close(
        next.frame.displacements()[3],
        STRESS / E + (STRESS - YIELD) / HARDENING + DT * RATE,
        2e-11,
    );
}
