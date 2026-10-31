//! G0/G1/G3: permanent creep changes native member extension, tension and
//! curvature. Manufactured coefficients are not measured alloy qualification.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::source::{NortonCreepSource, NortonSourceQuery, norton_requirements};
use fs_material::creep::uniaxial::UniaxialNorton;
use fs_material::creep::{CreepTolerances, NortonCreep, NortonCreepParameters};
use fs_material::elastic::IsotropicElastic;
use fs_material::state_point::MaterialPropertySelection;
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::creep::{FrameCreepError, FrameCreepTolerances};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel, SolidError};
use std::cell::Cell;

const E: f64 = 1e9;
const STRESS: f64 = 1e6;
const RATE: f64 = 1e-5;
const DT: f64 = 10.0;
const AREA: f64 = 5e-4;
const ARM: f64 = 0.05;
const FORCE: f64 = 2.0 * AREA * STRESS;

fn law(dt: f64, temperature: f64) -> UniaxialNorton {
    UniaxialNorton::new(
        NortonCreep::new(
            IsotropicElastic::new(E, 0.25, 0.03).unwrap(),
            NortonCreepParameters {
                reference_stress_pa: STRESS,
                reference_rate_s: RATE,
                stress_exponent: 1.0,
                activation_energy_j_mol: 20_000.0,
                reference_temperature_k: 300.0,
                temperature_range_k: [250.0, 500.0],
            },
            temperature,
            dt,
        )
        .unwrap(),
    )
}
fn section(law: &UniaxialNorton) -> SectionModel {
    SectionModel::Fiber(Section::new(
        [-ARM, ARM]
            .into_iter()
            .map(|y| Fiber {
                y,
                area: AREA,
                law: FiberLaw::Norton(law.clone()),
            })
            .collect(),
    ))
}
fn frame(law: &UniaxialNorton, axial_only: bool) -> Frame2d {
    let mut result = Frame2d::new(Geometry::Linear);
    result.add_node(0.0, 0.0);
    result.add_node(1.0, 0.0);
    result.fix(0, [true; 3]);
    if axial_only {
        result.fix(1, [false, true, true]);
    }
    result.add_element(0, 1, 3, 1.0, &|| section(law)).unwrap();
    result
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
fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={absolute:.3e}"
    );
}
fn axial_load(force: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, force, 0.0, 0.0]
}

#[test]
fn g1_constant_force_extends_the_native_member_and_leaves_permanent_set() {
    let mut accepted = frame(&law(DT, 300.0), true);
    for index in 1..=6 {
        let before = accepted.displacements().to_vec();
        let next = accepted
            .solve_creep_step(&axial_load(FORCE), controls(), || false)
            .unwrap();
        assert_eq!(accepted.displacements(), before);
        close(
            next.frame.displacements()[3],
            STRESS / E + index as f64 * DT * RATE,
            2e-11,
        );
        close(next.frame.element_forces(0)[0], FORCE, 1e-6);
        close(
            next.energy.physical_dissipation_j,
            2.0 * AREA * STRESS * DT * RATE,
            1e-8,
        );
        assert!(next.energy.balance_residual_j.abs() < controls().energy_j);
        assert_eq!(next.fibers.len(), 6);
        close(
            next.fibers.iter().map(|p| p.reference_volume_m3).sum(),
            2.0 * AREA,
            1e-16,
        );
        for p in 0..3 {
            for f in 0..2 {
                let (_, state) = next.frame.creep_fiber(0, p, f).unwrap();
                close(state.creep_strain(), index as f64 * DT * RATE, 2e-11);
                close(
                    state.equivalent_creep_strain(),
                    index as f64 * DT * RATE,
                    2e-11,
                );
            }
        }
        accepted = next.frame;
    }
    let released = accepted
        .solve_creep_step(&axial_load(0.0), controls(), || false)
        .unwrap();
    close(released.frame.displacements()[3], 6.0 * DT * RATE, 2e-11);
    close(released.frame.element_forces(0)[0], 0.0, 1e-6);
    assert!(released.energy.physical_dissipation_j < 1e-8);
}

#[test]
fn g1_held_member_relaxes_tension_and_heats_without_support_work() {
    let mut held = frame(&law(DT, 300.0), true)
        .solve_creep_step(&axial_load(FORCE), controls(), || false)
        .unwrap()
        .frame;
    held.fix(1, [true; 3]);
    let fixed = held.displacements().to_vec();
    let ratio = 1.0 / (1.0 + E * DT * RATE / STRESS);
    for index in 1..=10 {
        let step = held
            .solve_creep_step(&axial_load(0.0), controls(), || false)
            .unwrap();
        assert_eq!(step.frame.displacements(), fixed);
        close(
            step.frame.element_forces(0)[0],
            FORCE * ratio.powi(index),
            1e-6,
        );
        assert!(step.energy.physical_dissipation_j > 0.0);
        assert!(step.energy.numerical_dissipation_j > 0.0);
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
}

#[test]
fn g1_sustained_end_moment_changes_native_curvature_and_tip_shape() {
    let moment = 25.0;
    let inertia = 2.0 * AREA * ARM * ARM;
    let elastic_curvature = moment / (E * inertia);
    let creep_curvature_per_step = DT * RATE / STRESS * moment / inertia;
    let mut cantilever = frame(&law(DT, 300.0), false);
    let load = [0.0, 0.0, 0.0, 0.0, 0.0, moment];
    for index in 1..=5 {
        let step = cantilever
            .solve_creep_step(&load, controls(), || false)
            .unwrap();
        let curvature = elastic_curvature + index as f64 * creep_curvature_per_step;
        close(step.frame.displacements()[4], 0.5 * curvature, 2e-10);
        close(step.frame.displacements()[5], curvature, 2e-10);
        for &(axial, bending) in step.frame.section_deformations(0) {
            close(axial, 0.0, 1e-12);
            close(bending, curvature, 2e-10);
        }
        close(
            step.energy.physical_dissipation_j,
            moment * creep_curvature_per_step,
            1e-8,
        );
        let (_, bottom) = step.frame.creep_fiber(0, 1, 0).unwrap();
        let (_, top) = step.frame.creep_fiber(0, 1, 1).unwrap();
        assert!(bottom.creep_strain() > 0.0 && top.creep_strain() < 0.0);
        cantilever = step.frame;
    }
}

fn card(multiplier: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let values = [E, 0.25, STRESS, RATE * multiplier, 1.0, 20_000.0, 300.0];
    for (requirement, value) in norton_requirements().unwrap().iter().zip(values) {
        let dims = requirement.dims();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(requirement.name(), dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured frame creep fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured creep solid".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn source(card: &MaterialCard) -> NortonCreepSource<'_> {
    NortonCreepSource::new(
        card,
        NortonSourceQuery {
            lower: &QueryPoint::new().with("T", 250.0).unwrap(),
            upper: &QueryPoint::new().with("T", 500.0).unwrap(),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.03,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn g3_original_material_temperature_refresh_changes_the_continuing_member() {
    let material = card(1.0);
    let family = source(&material);
    let cold = UniaxialNorton::from_source(family.at_temperature(300.0, DT).unwrap());
    let accepted = frame(&cold, true)
        .solve_creep_step(&axial_load(FORCE), controls(), || false)
        .unwrap()
        .frame;
    let warm = UniaxialNorton::from_source(family.at_temperature(330.0, DT / 2.0).unwrap());
    let prepared = accepted
        .with_creep_laws(&vec![warm; 6], "T", || false)
        .unwrap();
    assert_eq!(prepared.displacements(), accepted.displacements());
    assert_eq!(
        prepared.creep_fiber(0, 0, 0).unwrap().1,
        accepted.creep_fiber(0, 0, 0).unwrap().1
    );
    let step = prepared
        .solve_creep_step(&axial_load(FORCE), controls(), || false)
        .unwrap();
    let thermal_factor = (20_000.0_f64 / 8.314_462_618_153_24 * (1.0 / 300.0 - 1.0 / 330.0)).exp();
    close(step.duration_s, DT / 2.0, 0.0);
    close(
        step.frame.displacements()[3] - accepted.displacements()[3],
        0.5 * DT * RATE * thermal_factor,
        2e-11,
    );
    let (next_law, _) = step.frame.creep_fiber(0, 2, 1).unwrap();
    let old = cold.material().unwrap();
    let next = next_law.material().unwrap();
    assert_eq!(old.card_identity(), next.card_identity());
    assert_eq!(next.state_coordinate("T"), Some(330.0));
    for (a, b) in old.properties().iter().zip(next.properties()) {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    let foreign = card(2.0);
    let foreign_law =
        UniaxialNorton::from_source(source(&foreign).at_temperature(330.0, DT).unwrap());
    assert!(
        accepted
            .with_creep_laws(&vec![foreign_law; 6], "T", || false)
            .is_err()
    );
    assert!(
        accepted
            .with_creep_laws(&vec![cold; 6], "pressure", || false)
            .is_err()
    );
}

#[test]
fn g0_refused_or_cancelled_physical_steps_preserve_all_accepted_history() {
    let accepted = frame(&law(DT, 300.0), true);
    let zero = accepted.displacements().to_vec();
    let history = accepted.creep_fiber(0, 0, 0).unwrap().1.clone();
    let mut too_small = controls();
    too_small.maximum_fiber_points = 5;
    assert!(
        accepted
            .solve_creep_step(&axial_load(FORCE), too_small, || false)
            .is_err()
    );
    let calls = Cell::new(0);
    let cancelled = accepted.solve_creep_step(&axial_load(FORCE), controls(), || {
        calls.set(calls.get() + 1);
        calls.get() == 12
    });
    assert!(matches!(cancelled, Err(FrameCreepError::Cancelled)));
    let mut mixed = vec![law(DT, 300.0); 6];
    mixed[3] = law(DT / 2.0, 300.0);
    let prepared = accepted.with_creep_laws(&mixed, "T", || false).unwrap();
    assert!(
        prepared
            .solve_creep_step(&axial_load(FORCE), controls(), || false)
            .is_err()
    );
    assert!(
        accepted
            .solve_creep_step(&axial_load(1000.0 * FORCE), controls(), || false)
            .is_err()
    );
    assert_eq!(accepted.displacements(), zero);
    assert_eq!(accepted.creep_fiber(0, 0, 0).unwrap().1, &history);
    let a = accepted
        .solve_creep_step(&axial_load(FORCE), controls(), || false)
        .unwrap();
    let b = accepted
        .solve_creep_step(&axial_load(FORCE), controls(), || false)
        .unwrap();
    assert_eq!(a.frame.displacements(), b.frame.displacements());
    assert_eq!(
        a.frame.creep_fiber(0, 0, 0).unwrap().1,
        b.frame.creep_fiber(0, 0, 0).unwrap().1
    );
}

#[test]
fn g0_untimed_or_inertial_routes_cannot_silently_advance_creep() {
    let mut accepted = frame(&law(DT, 300.0), true);
    let before = accepted.displacements().to_vec();
    let history = accepted.creep_fiber(0, 0, 0).unwrap().1.clone();
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
    assert_eq!(accepted.displacements(), before);
    assert_eq!(accepted.creep_fiber(0, 0, 0).unwrap().1, &history);
    let next = accepted
        .solve_creep_step(&axial_load(FORCE), controls(), || false)
        .unwrap();
    close(next.frame.displacements()[3], STRESS / E + DT * RATE, 2e-11);
}
