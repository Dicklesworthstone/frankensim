//! G0/G1/G3: directional sheet flow, springback and original Hill energy.
//! Coefficients are synthetic and do not calibrate a rolled alloy.

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::source::{
    HILL_COEFFICIENT_PROPERTIES, HILL_HARDENING_MODULUS_PROPERTY,
    HILL_REFERENCE_YIELD_STRESS_PROPERTY, ResolvedHill48, resolve_hill48,
};
use fs_material::plastic::hill::{Hill48, Hill48Coefficients, HillState};
use fs_material::plastic::section::energy::{NativeSectionEnergyLaw, SectionEnergyError};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_material::tensor::{contract, rotate, rotation};
use fs_qty::{Dims, Pressure};

const E: f64 = 1200.0;
const NU: f64 = 0.2;
const Y: f64 = 3.0;
const H: f64 = 40.0;

fn coefficients() -> Hill48Coefficients {
    Hill48Coefficients {
        f: 0.25,
        g: 0.9,
        h: 0.45,
        l: 1.4,
        m: 2.1,
        n: 1.8,
    }
}

fn native() -> Hill48 {
    Hill48::new(elastic(E, NU), Y, H, coefficients()).unwrap()
}

fn gates() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-9,
        yield_pa: 1e-8,
    }
}

#[test]
fn g1_directional_membrane_coupons_keep_local_flow_under_in_plane_rotation() {
    let c = coefficients();
    let alpha = 0.003;
    for axis in 0..2 {
        let flow = if axis == 0 {
            [c.g + c.h, -c.h, -c.g]
        } else {
            [-c.h, c.f + c.h, -c.f]
        };
        let k = flow[axis].sqrt();
        for sign in [-1.0, 1.0] {
            let sigma = sign * (Y + H * alpha) / k;
            let plastic: [f64; 6] = core::array::from_fn(|i| {
                if i < 3 {
                    sign * alpha * flow[i] / k
                } else {
                    0.0
                }
            });
            let local: [f64; 6] = core::array::from_fn(|i| {
                plastic[i]
                    + if i == axis {
                        sigma / E
                    } else if i < 3 {
                        -NU * sigma / E
                    } else {
                        0.0
                    }
            });
            let local_stress = core::array::from_fn(|i| if i == axis { sigma } else { 0.0 });
            for angle in [0.0, 0.47] {
                let frame = rotation([0.0, 0.0, 1.0], angle);
                let strain = rotate(&local, &frame);
                let stress = rotate(&local_stress, &frame);
                let law = OrientedHill48::new(native(), frame).unwrap();
                assert_eq!(law.physical_duration_s(), None);
                let section = section(law, 2);
                let old = section.initial_state();
                let next = section
                    .evaluate([strain[0], strain[1], 2.0 * strain[3], 0.0, 0.0, 0.0], &old)
                    .unwrap();
                let energy = section
                    .audit_energy(&old, &next.state, gates(), 1e-10)
                    .unwrap();
                for point in &next.state.points {
                    near(point.history.alpha, alpha, 2e-13);
                    for i in 0..6 {
                        near(point.history.plastic_strain[i], plastic[i], 2e-13);
                        near(point.strain[i], strain[i], 2e-13);
                    }
                }
                for (i, component) in [0, 1, 3].into_iter().enumerate() {
                    near(next.resultants[i], 0.2 * stress[component], 2e-11);
                    near(next.resultants[i + 3], 0.0, 2e-13);
                }
                near(energy.physical_dissipation_j_m2, 0.2 * Y * alpha, 1e-12);
                near(
                    energy.final_state.elastic_j_m2,
                    0.2 * sigma * sigma / (2.0 * E),
                    1e-12,
                );
                near(
                    energy.final_state.hardening_j_m2,
                    0.2 * 0.5 * H * alpha * alpha,
                    1e-12,
                );
                assert_eq!(energy.final_state.kinematic_hardening_j_m2, 0.0);
                near(energy.balance_residual_j_m2, 0.0, 1e-10);
                assert_eq!(old, section.initial_state());
            }
        }
    }
}

#[test]
fn g1_hill_partial_thickness_yield_retains_hardening_and_permanent_curvature() {
    let kinematic = 160.0;
    let law = native().with_linear_kinematic_hardening(kinematic).unwrap();
    assert_eq!(law.physical_duration_s(), None);
    let section = section(law, 8);
    let old = section.initial_state();
    let peak = [0.0, 0.0, 0.0, 0.04, 0.0, 0.0];
    let loaded = section.evaluate(peak, &old).unwrap();
    assert!(loaded.state.points.iter().any(|p| p.history.alpha == 0.0));
    assert!(loaded.state.points.iter().any(|p| p.history.alpha > 1e-4));
    let account = section
        .audit_energy(&old, &loaded.state, gates(), 1e-10)
        .unwrap();
    let mut heat = 0.0;
    let mut isotropic_storage = 0.0;
    let mut kinematic_storage = 0.0;
    for (point, state) in section.points().iter().zip(&loaded.state.points) {
        let history = &state.history;
        heat += point.weight_m * Y * history.alpha;
        isotropic_storage += point.weight_m * 0.5 * H * history.alpha * history.alpha;
        kinematic_storage += point.weight_m
            * 0.5
            * kinematic
            * contract(&history.plastic_strain, &history.plastic_strain);
    }
    assert!(heat > 0.0 && isotropic_storage > 0.0 && kinematic_storage > 0.0);
    near(account.physical_dissipation_j_m2, heat, 1e-12);
    near(account.final_state.hardening_j_m2, isotropic_storage, 1e-12);
    near(
        account.final_state.kinematic_hardening_j_m2,
        kinematic_storage,
        1e-12,
    );
    // Independent inverse of isotropic elastic bending stiffness gives the
    // stress-free springback; permanent anisotropic flow remains at each depth.
    let ei = E * 0.2_f64.powi(3) / 12.0;
    let mut free = peak;
    free[3] -= (loaded.resultants[3] - NU * loaded.resultants[4]) / ei;
    free[4] -= (loaded.resultants[4] - NU * loaded.resultants[3]) / ei;
    free[5] -= loaded.resultants[5] / (ei / (2.0 * (1.0 + NU)));
    assert!(free[3] > 0.0 && free[3] < peak[3]);
    let unloaded = section.evaluate(free, &loaded.state).unwrap();
    let release = section
        .audit_energy(&loaded.state, &unloaded.state, gates(), 1e-10)
        .unwrap();
    for resultant in unloaded.resultants {
        near(resultant, 0.0, 2e-11);
    }
    for (before, after) in loaded.state.points.iter().zip(&unloaded.state.points) {
        assert_eq!(before.history, after.history);
    }
    assert_eq!(release.physical_dissipation_j_m2, 0.0);
    assert_eq!(
        release.initial.hardening_j_m2,
        release.final_state.hardening_j_m2
    );
    assert_eq!(
        release.initial.kinematic_hardening_j_m2,
        release.final_state.kinematic_hardening_j_m2,
    );
    near(release.balance_residual_j_m2, 0.0, 1e-10);
}

fn sourced(frame: [[f64; 3]; 3]) -> ResolvedHill48 {
    let c = coefficients();
    let mut properties = vec![
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), NU),
        (HILL_REFERENCE_YIELD_STRESS_PROPERTY, Pressure::DIMS, Y),
        (HILL_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H),
    ];
    properties.extend(
        HILL_COEFFICIENT_PROPERTIES
            .into_iter()
            .zip([c.f, c.g, c.h, c.l, c.m, c.n])
            .map(|(name, value)| (name, Dims([0; 6]), value)),
    );
    let mut claims = ClaimSet::new();
    for (name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic layered Hill coupon".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored directional metal".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_hill48(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.05,
        frame,
    )
    .unwrap()
}

fn refuse_forged_history<L: NativeSectionEnergyLaw<State = HillState>>(law: L) {
    let section = section(law, 2);
    let old = section.initial_state();
    let next = section
        .evaluate([0.008, -0.001, 0.003, 0.01, 0.004, 0.006], &old)
        .unwrap();
    section
        .audit_energy(&old, &next.state, gates(), 1e-10)
        .unwrap();
    let saved = next.state.clone();
    let mut incoming = saved.clone();
    for point in &mut incoming.points {
        point.history = HillState::default();
    }
    // Same represented strain with erased plastic memory is overstressed at
    // the beginning. It must refuse instead of being repaired by a new return.
    assert!(matches!(
        section.audit_energy(&incoming, &next.state, gates(), 1e-10),
        Err(SectionEnergyError::Native { .. })
    ));
    let mut substituted = saved.clone();
    substituted.points[0].history.alpha += 1e-4;
    assert!(matches!(
        section.audit_energy(&old, &substituted, gates(), 1e-10),
        Err(SectionEnergyError::State { .. })
    ));
    assert_eq!(next.state, saved);
}

#[test]
fn g3_original_hill_source_preserves_section_response_history_and_energy() {
    let frame = rotation([0.0, 0.0, 1.0], 0.37);
    let resolved = sourced(frame);
    let saved = resolved.clone();
    assert_eq!(resolved.physical_duration_s(), None);
    let original = section(OrientedHill48::new(native(), frame).unwrap(), 3);
    let sourced = section(resolved.clone(), 3);
    let mut old = original.initial_state();
    for generalized in [
        [0.007, -0.001, 0.003, 0.025, -0.01, 0.013],
        [-0.003, 0.002, -0.001, -0.01, 0.005, -0.007],
    ] {
        let a = original.evaluate(generalized, &old).unwrap();
        let b = sourced.evaluate(generalized, &old).unwrap();
        assert_eq!(a.state, b.state);
        assert_eq!(a.resultants, b.resultants);
        assert_eq!(a.tangent, b.tangent);
        assert_eq!(
            original
                .audit_energy(&old, &a.state, gates(), 1e-10)
                .unwrap(),
            sourced
                .audit_energy(&old, &b.state, gates(), 1e-10)
                .unwrap(),
        );
        old = a.state;
    }
    assert_eq!(resolved, saved);
    refuse_forged_history(resolved);
}

#[test]
fn g0_hill_energy_refuses_inadmissible_old_or_substituted_new_history() {
    let law = native().with_linear_kinematic_hardening(160.0).unwrap();
    refuse_forged_history(law.clone());
    refuse_forged_history(OrientedHill48::new(law, rotation([0.0, 0.0, 1.0], 0.42)).unwrap());
}
