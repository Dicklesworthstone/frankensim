//! Source-selected fast Maxwell memory in the actual acoustic wall owner.
use super::super::{
    AxialConstraint, LocalMaterialWall, MaterialWallAssignment, WallMaterialModel,
    source::WallMaterialLossSource,
};
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};
use fs_qty::Density;
use fs_vfit::waveguide::network::NetworkSegment;

const BAND: [f64; 2] = [20.0, 4000.0];
const DT: f64 = 1.0 / 48_000.0;

fn selection(exact: bool) -> WallMaterialLossSource {
    let mode = if exact { "maxwell-exact" } else { "maxwell" };
    WallMaterialLossSource::read(
        &mut format!("{mode} 2 1.0 0.8 arm_e1 arm_tau1 arm_e2 arm_tau2").split_ascii_whitespace(),
    )
    .unwrap()
}

fn material() -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, 1e6),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("arm_e1", Pressure::DIMS, 8e5),
        ("arm_tau1", Dims([0, 0, 1, 0, 0, 0]), 1e-6),
        ("arm_e2", Pressure::DIMS, 4e5),
        ("arm_tau2", Dims([0, 0, 1, 0, 0, 0]), 2e-7),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic fast-relaxing wall; not measured instrument material".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-fast-wall".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &selection(true).requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn bound(material: &ResolvedMaterialStatePoint, exact: bool) -> BoundMaterialWall {
    let source = MaterialWallAssignment {
        node: 1,
        material: material.clone(),
        model: WallMaterialModel::Isotropic,
        thickness_m: 0.0002,
        axial: AxialConstraint::StressFree,
        loss: selection(exact).resolve(material, BAND).unwrap(),
    };
    let local = LocalMaterialWall::new(
        source.material.clone(),
        source.thickness_m,
        source.axial,
        source.loss.clone(),
        BAND,
    )
    .unwrap();
    let (radius, length) = (0.006, 0.15);
    let patch = local.base_patch(radius, length).unwrap();
    let load = local.load(patch).unwrap();
    BoundMaterialWall {
        source,
        radius_m: radius,
        length_m: length,
        band_hz: BAND,
        patch,
        load,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * scale.max(1e-30),
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

#[test]
fn exact_source_selection_admits_fast_memory_and_retains_the_membrane_clock() {
    let material = material();
    let midpoint = bound(&material, false);
    let exact = bound(&material, true);
    assert_eq!(
        selection(false).requirements(),
        selection(true).requirements()
    );
    assert_eq!(exact.patch(), midpoint.patch());
    assert_eq!(exact.load().base(), midpoint.load().base());
    assert_eq!(
        exact.load().compliance_terms(),
        midpoint.load().compliance_terms()
    );
    assert_eq!(
        midpoint.load().sampling(),
        RelaxationSampling::BilinearMidpoint
    );
    assert_eq!(exact.load().sampling(), RelaxationSampling::ExactHeldFlow);
    assert!(midpoint.validate_time_step(DT).is_err());
    exact.validate_time_step(DT).unwrap();
    // Exact stress memory does not make the midpoint inertial membrane exact.
    assert!(exact.validate_time_step(0.001).is_err());
    let WallMaterialLoss::Maxwell(source) = &exact.source.loss else {
        unreachable!()
    };
    assert_eq!(source.source_identity(), Some(material.identity()));
    // The separate distributed FE owner still uses midpoint stress memory.
    assert!(
        source
            .validate_local_time_step(exact.patch.wall, DT)
            .is_err()
    );
    for dt in [0.0, -DT, f64::NAN, f64::INFINITY] {
        assert!(exact.validate_time_step(dt).is_err());
    }
}

#[test]
fn exact_source_sampling_survives_pinned_refresh_and_cannot_change_live_history() {
    let material = material();
    let wall = bound(&material, true);
    let WallMaterialLoss::Maxwell(source) = &wall.source.loss else {
        unreachable!()
    };
    let refreshed = source.at_material_pinned(&material, &material).unwrap();
    assert_eq!(source, &refreshed);
    assert_eq!(
        source.history_correspondence(&refreshed).unwrap(),
        [(0, 0), (1, 1)]
    );
    let changed = refreshed.with_network_sampling(RelaxationSampling::BilinearMidpoint);
    assert!(source.history_correspondence(&changed).is_err());
}

#[test]
fn fast_sourced_wall_memory_closes_motion_heat_and_the_full_network_work() {
    let wall = bound(&material(), true);
    wall.validate_time_step(DT).unwrap();
    let mut network = WaveguideNetwork::new(
        &[
            NetworkNode::Inlet,
            wall.shunt(),
            NetworkNode::Termination { reflection: -0.6 },
        ],
        &[
            NetworkSegment {
                nodes: [0, 1],
                one_way_samples: 7,
                impedance_pa_s_m3: 3.5e6,
            },
            NetworkSegment {
                nodes: [1, 2],
                one_way_samples: 11,
                impedance_pa_s_m3: 3.5e6,
            },
        ],
        DT,
        1 << 20,
    )
    .unwrap();
    let mut old = wall.observe_waveguide(&network).unwrap();
    let mut total_heat = 0.0;
    for sample in 0..900 {
        let old_pressures = network.load_compliance_pressures(1).unwrap().to_vec();
        let old_energy = network.stored_energy_j();
        let frame = network.step(if sample < 60 { 10.0 } else { 0.0 }).unwrap();
        let new = wall.observe_waveguide(&network).unwrap();
        let node = network.node_frame(1).unwrap();
        let pressures = network.load_compliance_pressures(1).unwrap();
        let (mut mean_pressure, mut heat) = (0.0, 0.0);
        for (j, term) in wall.load().compliance_terms().iter().enumerate() {
            let rate = term.rate_per_s;
            let resistance = term.resistance_pa_s_m3;
            let old_p = old_pressures[j];
            let target = resistance * node.load_flow_m3_s;
            let difference = old_p - target;
            let decay = (-rate * DT).exp();
            let decrement = -(-rate * DT).exp_m1();
            let expected = target + decay * difference;
            close(pressures[j], expected, old_p.abs() + target.abs());
            mean_pressure += target + difference * decrement / (rate * DT);
            // Independent analytic integral of p(t)^2/R. The fixture's large
            // rate*DT keeps this direct formula well conditioned.
            heat += (target * target * DT
                + 2.0 * target * difference * decrement / rate
                + difference * difference * -(-2.0 * rate * DT).exp_m1() / (2.0 * rate))
                / resistance;
        }
        assert!(heat >= 0.0);
        let inertia = wall.patch.wall.surface_density * (new.velocity_m_s - old.velocity_m_s) / DT;
        let spring =
            wall.patch.wall.stiffness_per_area * 0.5 * (new.displacement_m + old.displacement_m);
        close(
            inertia + spring + mean_pressure,
            node.pressure_pa,
            inertia.abs() + spring.abs() + mean_pressure.abs() + node.pressure_pa.abs(),
        );
        close(
            wall.patch.area_m2 * 0.5 * (new.velocity_m_s + old.velocity_m_s),
            node.load_flow_m3_s,
            node.load_flow_m3_s.abs(),
        );
        close(heat, node.absorbed_energy_j, heat);
        let stored = new.stored_energy_j + wall.relaxation_stored_energy_j(&network).unwrap();
        close(stored, node.stored_energy_j, stored);
        let work = node.pressure_pa * node.load_flow_m3_s * DT;
        close(
            node.storage_change_j + heat,
            work,
            node.storage_change_j.abs() + heat + work.abs(),
        );
        close(
            frame.balance_residual_j(),
            0.0,
            old_energy + frame.stored_energy_j + frame.inlet_work_j.abs(),
        );
        total_heat += heat;
        old = new;
    }
    assert!(total_heat > 0.0);
}
