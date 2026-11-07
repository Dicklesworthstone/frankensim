//! Original ordered-interface spectra reach the native coupled load line.
use fs_couple::cohesive_contact::{CohesiveContactConfig, CohesiveLoadLine};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::relaxation::{ViscoelasticBondLimits, source::*};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, Length, Pressure, Time};

fn limits() -> ViscoelasticBondLimits {
    ViscoelasticBondLimits {
        maximum_jump_m: 1e-4,
        maximum_rate_m_s: 1.0,
        maximum_elapsed_s: 2.0,
    }
}
fn card(missing_shear: bool, wrong_units: bool) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    let mut parameters = vec![(BONDLINE_THICKNESS_PROPERTY.to_owned(), Length::DIMS, 0.001)];
    for (prefix, e0, e1, tau) in [
        (NORMAL_RELAXATION_PREFIX, 1e6, 2e6, 0.01),
        (SHEAR_RELAXATION_PREFIX, 4e5, 8e5, 0.025),
    ] {
        if missing_shear && prefix == SHEAR_RELAXATION_PREFIX {
            continue;
        }
        parameters.extend([
            (format!("{prefix}_branch_count"), Dims::NONE, 1.0),
            (format!("{prefix}_equilibrium_modulus"), Pressure::DIMS, e0),
            (format!("{prefix}_modulus_1"), Pressure::DIMS, e1),
            (
                format!("{prefix}_time_1"),
                if wrong_units {
                    Length::DIMS
                } else {
                    Time::DIMS
                },
                tau,
            ),
        ]);
    }
    for (name, dims, value) in parameters {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured intact joint relaxation; not measured adhesive data"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "specified prepared substrate".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "fixed intact thin adhesive".into(),
            third_body: None,
            environment: "fixed temperature".into(),
            history: "declared cure and normal restraint".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
#[test]
fn sourced_relaxation_drives_native_bonded_motion_and_retains_domain_refusals() {
    let card = card(false, false);
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    let source = resolve_viscoelastic_bond(
        &card,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        limits(),
    )
    .unwrap();
    assert_eq!(source.resolved().card_identity(), card.content_hash());
    assert_eq!(source.resolved().properties().len(), 9);
    let pins = source
        .resolved()
        .properties()
        .iter()
        .map(|p| {
            (
                p.requirement().name().to_owned(),
                p.answer().receipt.selected,
            )
        })
        .collect::<Vec<_>>();
    let pinned = resolve_viscoelastic_bond(
        &card,
        &query,
        MaterialPropertySelection::PinnedByProperty(pins.clone()),
        limits(),
    )
    .unwrap();
    assert_eq!(pinned.law(), source.law());
    for (a, b) in pinned
        .resolved()
        .properties()
        .iter()
        .zip(source.resolved().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    let mut foreign_pins = pins;
    foreign_pins.push(("unowned_relaxation_coefficient".into(), foreign_pins[0].1));
    assert!(
        resolve_viscoelastic_bond(
            &card,
            &query,
            MaterialPropertySelection::PinnedByProperty(foreign_pins),
            limits()
        )
        .is_err()
    );
    assert_eq!(source.surfaces()[0], *card.surface_a());
    let dt = 0.001;
    let area = 1e-4;
    let compliance = 0.5 * dt * dt;
    let config = CohesiveContactConfig {
        max_iterations: 8,
        max_material_evaluations: 16,
        max_backtracks: 8,
        jump_absolute_tolerance_m: 1e-14,
        jump_relative_tolerance: 1e-11,
        maximum_jump_increment_m: 1e-4,
        maximum_force_n: 100.0,
        material_energy_tolerance_j_m2: 1e-9,
    };
    let mut solver = CohesiveLoadLine::new(
        area,
        [
            [compliance, 0.0, 0.0],
            [0.0, compliance, 0.0],
            [0.0, 0.0, compliance],
        ],
        config,
    )
    .unwrap();
    let prepared = source.law().prepare(dt).unwrap();
    let mut point = source.law().bonded_point();
    let mut velocity = [0.0; 3];
    let load = [0.5, -0.2, 0.1];
    let mut work = 0.0;
    let mut heat = 0.0;
    for _ in 0..1000 {
        let old = point.jump_m();
        let free = core::array::from_fn(|i| old[i] + dt * velocity[i] + compliance * load[i]);
        let step = solver
            .trial_relaxation(&point, &prepared, free, || true)
            .unwrap();
        let new = step.material.point.jump_m();
        for i in 0..3 {
            work += load[i] * (new[i] - old[i]);
            velocity[i] = 2.0 * (new[i] - old[i]) / dt - velocity[i];
        }
        heat += area * step.material.viscous_dissipated_j_m2;
        point = step.material.point;
    }
    assert!(heat > 0.0);
    let kinetic = 0.5 * velocity.iter().map(|v| v * v).sum::<f64>();
    assert!((work - kinetic - area * point.stored_energy_j_m2() - heat).abs() < 1e-10);
    assert!((point.jump_m()[0] - 5e-6).abs() < 1e-9);
    assert!(
        point
            .trial(&source.law().prepare(1.1).unwrap(), point.jump_m(), 1e-9)
            .is_err()
    );
    assert!(
        resolve_viscoelastic_bond(
            &card,
            &QueryPoint::new().with("T", 321.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            limits()
        )
        .is_err()
    );
}
#[test]
fn incomplete_and_dimensionally_wrong_interface_spectra_do_not_become_bulk_defaults() {
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    for interface in [card(true, false), card(false, true)] {
        assert!(
            resolve_viscoelastic_bond(
                &interface,
                &query,
                MaterialPropertySelection::SingleClaimOnly,
                limits()
            )
            .is_err()
        );
    }
}
