//! Bead sj31i.27 battery: the wedge fixture's π semantics are honest.
//! The v1 fixture named kinematic dims "viscosity" while carrying air's
//! DYNAMIC viscosity value, silently building V·L/ν with no density.
//! The v2 basis carries role-checked density and dynamic viscosity, so
//! Buckingham derives the true Reynolds group ρVL/μ; role-incompatible
//! substitutions refuse typed; the corrupted v1 basis diverges
//! structurally instead of aliasing.
#![cfg(feature = "flywheel-e2e")]

use std::collections::BTreeMap;

use fs_flywheel_e2e::{
    KINEMATIC_VISCOSITY_DIMS, WEDGE_ROLE_DENSITY, WEDGE_ROLE_DYNAMIC_VISCOSITY, WEDGE_ROLE_LENGTH,
    WEDGE_ROLE_VELOCITY, WedgeMaterial, WedgeRoleError, wedge_descriptor,
};
use fs_ledger::tombstone::{Descriptor, pi_distance};
use fs_qty::{Dims, QtyAny};

/// The two textbook Reynolds forms agree through the CHECKED quantity
/// algebra: V·L/ν (ν derived as μ/ρ) and ρ·V·L/μ both cancel to
/// dimensionless and match within a certified rounding bound.
#[test]
fn equivalent_reynolds_forms_agree_within_certified_bound() {
    for material in [WedgeMaterial::air(), WedgeMaterial::water()] {
        let velocity = QtyAny::new(20.0, WEDGE_ROLE_VELOCITY.1);
        let length = QtyAny::new(0.1, WEDGE_ROLE_LENGTH.1);

        let nu = material
            .kinematic_viscosity()
            .expect("role-checked material derives kinematic viscosity");
        assert_eq!(nu.dims, KINEMATIC_VISCOSITY_DIMS);

        let re_kinematic = velocity
            .try_mul(length)
            .expect("V*L stays in range")
            .try_div(nu)
            .expect("V*L/nu stays in range");
        let re_dynamic = material
            .density
            .try_mul(velocity)
            .expect("rho*V stays in range")
            .try_mul(length)
            .expect("rho*V*L stays in range")
            .try_div(material.dynamic_viscosity)
            .expect("rho*V*L/mu stays in range");

        assert_eq!(re_kinematic.dims, Dims::NONE, "V*L/nu is dimensionless");
        assert_eq!(re_dynamic.dims, Dims::NONE, "rho*V*L/mu is dimensionless");
        let relative = ((re_kinematic.value - re_dynamic.value) / re_dynamic.value).abs();
        assert!(
            relative <= 4.0 * f64::EPSILON,
            "{}: Re forms diverge beyond the rounding bound: {} vs {} (rel {relative:e})",
            material.label,
            re_kinematic.value,
            re_dynamic.value
        );
    }
}

/// The descriptor's π basis IS the executable dimensional derivation:
/// four role-named parameters over three base dimensions yield exactly
/// one dimensionless group whose value is the Reynolds number (up to
/// the Buckingham sign convention), and replay is deterministic.
#[test]
fn wedge_pi_group_is_the_reynolds_number_and_replays() {
    let air = WedgeMaterial::air();
    let descriptor = wedge_descriptor("cht-wedge bracket", 20.0, 0.1, &air);
    let signature = descriptor.pi_signature().expect("v2 basis admits");

    assert_eq!(
        signature.basis,
        vec![
            "density".to_string(),
            "dynamic_viscosity".to_string(),
            "length".to_string(),
            "velocity".to_string(),
        ],
        "the basis binds every semantic role, density included"
    );
    assert_eq!(
        signature.groups.len(),
        1,
        "four parameters over three base dimensions give exactly one group"
    );

    let re: f64 = 1.225 * 20.0 * 0.1 / 1.8e-5;
    let group_log10 = signature.groups[0].1;
    assert!(
        (group_log10.abs() - re.log10()).abs() < 1e-9,
        "the single group is Re (or its reciprocal): |log10| {} vs {}",
        group_log10.abs(),
        re.log10()
    );

    let replay = wedge_descriptor("cht-wedge bracket", 20.0, 0.1, &air)
        .pi_signature()
        .expect("replay admits");
    assert_eq!(signature, replay, "π signatures replay bitwise");
}

/// Metamorphic rescaling: scaling ρ and μ by the same factor leaves the
/// physics (ν and Re) invariant — bitwise for a power-of-two factor.
#[test]
fn common_density_viscosity_rescaling_leaves_the_group_invariant() {
    let air = WedgeMaterial::air();
    let scale = 1024.0; // exact power of two: quotients are bitwise-stable
    let scaled = WedgeMaterial::new(
        "air-rescaled",
        QtyAny::new(air.density.value * scale, WEDGE_ROLE_DENSITY.1),
        QtyAny::new(
            air.dynamic_viscosity.value * scale,
            WEDGE_ROLE_DYNAMIC_VISCOSITY.1,
        ),
    )
    .expect("rescaled properties keep their roles");

    let nu = air.kinematic_viscosity().expect("air nu");
    let nu_scaled = scaled.kinematic_viscosity().expect("scaled nu");
    assert_eq!(
        nu.value.to_bits(),
        nu_scaled.value.to_bits(),
        "nu = mu/rho is invariant under a common power-of-two rescale"
    );

    let base = wedge_descriptor("w", 20.0, 0.1, &air)
        .pi_signature()
        .expect("base");
    let rescaled = wedge_descriptor("w", 20.0, 0.1, &scaled)
        .pi_signature()
        .expect("rescaled");
    assert_eq!(base.groups[0].0, rescaled.groups[0].0, "same exponents");
    assert!(
        (base.groups[0].1 - rescaled.groups[0].1).abs() < 1e-12,
        "the Reynolds group value is rescale-invariant"
    );
}

/// Air and water are the same physics family (identical basis and
/// exponents) at genuinely different group values.
#[test]
fn air_and_water_share_structure_but_not_value() {
    let air = wedge_descriptor("w", 20.0, 0.1, &WedgeMaterial::air())
        .pi_signature()
        .expect("air");
    let water = wedge_descriptor("w", 20.0, 0.1, &WedgeMaterial::water())
        .pi_signature()
        .expect("water");
    assert_eq!(air.basis, water.basis);
    assert_eq!(air.groups[0].0, water.groups[0].0);
    let distance = pi_distance(&air, &water).expect("same structure is comparable");
    assert!(
        distance > 0.5,
        "air and water Reynolds numbers differ by decades: {distance}"
    );
}

/// Role incompatibility refuses TYPED: the v1 bug's exact shape
/// (kinematic dims offered as dynamic viscosity), a density with wrong
/// dims, and every single-exponent mutation of either property role.
#[test]
fn role_incompatible_substitutions_refuse_typed() {
    // The v1 bug reconstructed: 1.8e-5 with KINEMATIC dims cannot enter
    // the dynamic-viscosity role, however plausible the magnitude.
    let v1_bug = WedgeMaterial::new(
        "v1-corrupted",
        QtyAny::new(1.225, WEDGE_ROLE_DENSITY.1),
        QtyAny::new(1.8e-5, KINEMATIC_VISCOSITY_DIMS),
    );
    assert_eq!(
        v1_bug,
        Err(WedgeRoleError {
            role: "dynamic_viscosity",
            expected: WEDGE_ROLE_DYNAMIC_VISCOSITY.1,
            found: KINEMATIC_VISCOSITY_DIMS,
        })
    );

    // A kinematic value offered as density refuses on the density role.
    let wrong_density = WedgeMaterial::new(
        "bad-density",
        QtyAny::new(1.0, KINEMATIC_VISCOSITY_DIMS),
        QtyAny::new(1.8e-5, WEDGE_ROLE_DYNAMIC_VISCOSITY.1),
    );
    assert!(matches!(
        wrong_density,
        Err(WedgeRoleError {
            role: "density",
            ..
        })
    ));

    // Every single-exponent mutation of each property role refuses.
    for (role_dims, other_dims, mutated_role) in [
        (
            WEDGE_ROLE_DENSITY.1,
            WEDGE_ROLE_DYNAMIC_VISCOSITY.1,
            "density",
        ),
        (
            WEDGE_ROLE_DYNAMIC_VISCOSITY.1,
            WEDGE_ROLE_DENSITY.1,
            "dynamic_viscosity",
        ),
    ] {
        for position in 0..6 {
            for delta in [-1i8, 1] {
                let mut mutated = role_dims;
                mutated.0[position] += delta;
                let (density, viscosity) = if mutated_role == "density" {
                    (QtyAny::new(1.0, mutated), QtyAny::new(1.0, other_dims))
                } else {
                    (QtyAny::new(1.0, other_dims), QtyAny::new(1.0, mutated))
                };
                let refused = WedgeMaterial::new("mutant", density, viscosity);
                assert!(
                    matches!(&refused, Err(e) if e.role == mutated_role && e.found == mutated),
                    "exponent {position} delta {delta} on {mutated_role} must refuse: {refused:?}"
                );
            }
        }
    }
}

/// The corrupted v1 basis diverges STRUCTURALLY from v2: different
/// parameter roles mean no π comparison at all (no collision claim in
/// either direction), so stale tombstone signatures cannot alias the
/// corrected physics.
#[test]
fn v1_corrupted_basis_diverges_instead_of_aliasing() {
    let mut v1_params = BTreeMap::new();
    v1_params.insert(
        "velocity".to_string(),
        QtyAny::new(20.0, WEDGE_ROLE_VELOCITY.1),
    );
    v1_params.insert("length".to_string(), QtyAny::new(0.1, WEDGE_ROLE_LENGTH.1));
    v1_params.insert(
        "viscosity".to_string(),
        QtyAny::new(1.8e-5, KINEMATIC_VISCOSITY_DIMS),
    );
    let v1 = Descriptor {
        name: "cht-wedge bracket".to_string(),
        params: v1_params,
    }
    .pi_signature()
    .expect("the v1 shape still forms a (corrupted) signature");

    let v2 = wedge_descriptor("cht-wedge bracket", 20.0, 0.1, &WedgeMaterial::air())
        .pi_signature()
        .expect("v2");

    assert_ne!(v1.basis, v2.basis, "the role rename IS the crosswalk");
    assert_eq!(
        pi_distance(&v1, &v2),
        None,
        "structurally different bases make no distance claim either way"
    );
}
