//! Synthetic source cards test constitutive propagation, not measured flute fidelity.
#[path = "material_duct_wall/flue.rs"]
mod flue;
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::wall::material::{
    AxialConstraint, MaterialWallAssignment, WallMaterialLoss, YOUNG_VISCOSITY_PROPERTY,
    bind_material_walls,
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, Dims, DynViscosity, Pressure, QuantitySpec};

const RATE: u32 = 48_000;
const BAND: [f64; 2] = [20.0, 4000.0];

fn state(
    properties: &[(&str, Dims, f64)],
    band: Option<[f64; 2]>,
    curve: Option<&str>,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(key, dims, value) in properties {
        let mut validity = ValidityDomain::unconstrained().with("T", 290.0, 300.0);
        if let Some([lo, hi]) = band {
            validity = validity.with("omega", lo, hi);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: if curve == Some(key) {
                    PropertyValue::Curve {
                        abscissa: "omega".into(),
                        abscissa_dims: Dims([0, 0, -1, 0, 0, 0]),
                        knots: vec![(1.0, value), (30_000.0, 2.0 * value)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if curve == Some(key) {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic isotropic cylinder; not measured material data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                key,
                QuantitySpec::dimensional(dims),
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-wall".into(),
            phase: "solid".into(),
            process: "manufactured test".into(),
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
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn properties(rho: f64, young: f64, eta: f64) -> [(&'static str, Dims, f64); 4] {
    [
        ("density", Density::DIMS, rho),
        ("young_modulus", Pressure::DIMS, young),
        ("poisson_ratio", Dims([0; 6]), 0.3),
        (YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS, eta),
    ]
}
fn assignment(rho: f64, young: f64, eta: f64) -> MaterialWallAssignment {
    MaterialWallAssignment {
        node: 1,
        material: state(&properties(rho, young, eta), None, None),
        thickness_m: 0.0002,
        axial: AxialConstraint::StressFree,
        loss: WallMaterialLoss::KelvinVoigt,
    }
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn graph() -> TubeNetworkSpec {
    TubeNetworkSpec {
        nodes: vec![
            NetworkNode::Inlet,
            NetworkNode::Junction,
            NetworkNode::Termination { reflection: -0.9 },
        ],
        sections: (0..2)
            .map(|i| TubeSection {
                nodes: [i, i + 1],
                length_m: 0.15,
                radius_m: 0.006,
                max_length_error_m: 0.005,
            })
            .collect(),
        sound_speed_m_s: gas().sound_speed,
        max_wave_memory_bytes: 1_000_000,
    }
}
fn bound(
    a: MaterialWallAssignment,
) -> fs_couple::bernoulli_aperture::wall::material::BoundMaterialWall {
    bind_material_walls(graph(), vec![a], BAND, 1, &CancelGate::new_clock_free())
        .unwrap()
        .1
        .remove(0)
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * (a.abs() + b.abs()).max(1e-30),
        "{a} != {b}"
    );
}

#[test]
fn original_cards_and_geometry_drive_each_membrane_coefficient() {
    let a = assignment(1000.0, 1e6, 2.0);
    let w = bound(a.clone());
    assert_eq!(w.source(), &a);
    assert_eq!(w.band_hz(), BAND);
    assert_eq!(w.radius_m(), 0.006);
    assert_eq!(w.length_m(), 0.15);
    close(w.patch().area_m2, core::f64::consts::TAU * 0.006 * 0.15);
    close(w.patch().wall.surface_density, 0.2);
    close(
        w.patch().wall.stiffness_per_area,
        1e6 * 0.0002 / 0.006_f64.powi(2),
    );
    close(w.patch().wall.resistance, 2.0 * 0.0002 / 0.006_f64.powi(2));
    let mass = bound(assignment(2000.0, 1e6, 2.0));
    close(
        mass.patch().wall.surface_density,
        2.0 * w.patch().wall.surface_density,
    );
    assert_eq!(
        mass.patch().wall.stiffness_per_area,
        w.patch().wall.stiffness_per_area
    );
    let stiffness = bound(assignment(1000.0, 4e6, 2.0));
    close(
        stiffness.patch().wall.stiffness_per_area,
        4.0 * w.patch().wall.stiffness_per_area,
    );
    assert_eq!(
        stiffness.patch().wall.surface_density,
        w.patch().wall.surface_density
    );
    let mut fixed = a.clone();
    fixed.axial = AxialConstraint::StrainFixed;
    let fixed = bound(fixed);
    close(
        fixed.patch().wall.stiffness_per_area,
        w.patch().wall.stiffness_per_area / 0.91,
    );
    close(
        fixed.patch().wall.resistance,
        w.patch().wall.resistance / 0.91,
    );
    // Elastic is an explicit loss choice, not missing-data fallback.
    let mut elastic = a;
    elastic.loss = WallMaterialLoss::Elastic;
    elastic.material = state(&properties(1000.0, 1e6, 2.0)[..3], None, None);
    assert_eq!(bound(elastic).patch().wall.resistance, 0.0);
}

#[test]
fn adjacent_assignments_preserve_addresses_and_do_not_double_count_area() {
    let mut g = graph();
    g.nodes.insert(2, NetworkNode::Junction);
    g.sections.push(TubeSection {
        nodes: [2, 3],
        ..g.sections[0]
    });
    let original = g.clone();
    let a = assignment(1000.0, 1e6, 2.0);
    let mut b = a.clone();
    b.node = 2;
    let (g, walls) =
        bind_material_walls(g, vec![a, b], BAND, 2, &CancelGate::new_clock_free()).unwrap();
    assert_eq!(g.sections, original.sections);
    assert_eq!(g.nodes.len(), original.nodes.len());
    assert_eq!(g.nodes[0], original.nodes[0]);
    assert_eq!(g.nodes[3], original.nodes[3]);
    assert_eq!(g.sound_speed_m_s, original.sound_speed_m_s);
    assert_eq!(g.max_wave_memory_bytes, original.max_wave_memory_bytes);
    // Three sections: two interior strips cover the middle plus half of each end.
    close(
        walls.iter().map(|w| w.patch().area_m2).sum(),
        core::f64::consts::TAU * 0.006 * 0.3,
    );
    for w in walls {
        assert_eq!(g.nodes[w.source().node], w.patch().shunt().unwrap());
    }
}

#[test]
fn incompatible_material_topology_domain_or_budget_never_silently_falls_back() {
    let gate = CancelGate::new_clock_free();
    let refuses = |g, a| assert!(bind_material_walls(g, vec![a], BAND, 1, &gate).is_err());
    let a = assignment(1000.0, 1e6, 2.0);
    let mut missing = a.clone();
    missing.material = state(&properties(1000.0, 1e6, 2.0)[..3], None, None);
    refuses(graph(), missing);
    let mut wrong = properties(1000.0, 1e6, 2.0);
    wrong[3].1 = Pressure::DIMS;
    let mut b = a.clone();
    b.material = state(&wrong, None, None);
    refuses(graph(), b);
    for (idx, value) in [(0, 0.0), (1, 0.0), (2, 0.5), (2, -1.0), (3, -1.0)] {
        let mut p = properties(1000.0, 1e6, 2.0);
        p[idx].2 = value;
        let mut b = a.clone();
        b.material = state(&p, None, None);
        refuses(graph(), b);
    }
    let mut b = a.clone();
    b.thickness_m = 0.001;
    refuses(graph(), b);
    let mut b = a.clone();
    b.material = state(&properties(1000.0, 1e6, 2.0), Some([1.0, 2000.0]), None);
    refuses(graph(), b);
    let mut covered = a.clone();
    covered.material = state(
        &properties(1000.0, 1e6, 2.0),
        Some(BAND.map(|f| core::f64::consts::TAU * f)),
        None,
    );
    assert_eq!(bound(covered).band_hz(), BAND);
    let mut b = a.clone();
    b.material = state(&properties(1000.0, 1e6, 2.0), None, Some("young_modulus"));
    refuses(graph(), b);
    let mut g = graph();
    g.sections[1].radius_m *= 2.0;
    refuses(g, a.clone());
    let mut g = graph();
    g.nodes[1] = bound(a.clone()).patch().shunt().unwrap();
    refuses(g, a.clone());
    let mut g = graph();
    g.sections.push(g.sections[1]);
    refuses(g, a.clone());
    assert!(bind_material_walls(graph(), vec![a.clone(), a.clone()], BAND, 2, &gate).is_err());
    assert!(bind_material_walls(graph(), vec![a.clone()], BAND, 0, &gate).is_err());
    assert!(bind_material_walls(graph(), vec![a.clone()], [0.0, 4000.0], 1, &gate).is_err());
    gate.request();
    refuses(graph(), a);
}

#[test]
fn material_wall_motion_closes_the_existing_pressure_flow_and_energy_equations() {
    let dt = 1.0 / f64::from(RATE);
    let mut signals = Vec::new();
    for young in [1e6, 4e6] {
        let (g, walls) = bind_material_walls(
            graph(),
            vec![assignment(1000.0, young, 2.0)],
            BAND,
            1,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let patch = walls[0].patch();
        let (mut n, _) = g
            .realize_waveguide(&[gas(), gas()], dt, &CancelGate::new_clock_free())
            .unwrap();
        let mut old = patch.observe_waveguide(&n, 1).unwrap();
        let mut trace = Vec::new();
        let mut peak = 0.0_f64;
        let mut total_loss = 0.0;
        for i in 0..1024 {
            let f = n.step(if i < 20 { 1.0 } else { 0.0 }).unwrap();
            let new = patch.observe_waveguide(&n, 1).unwrap();
            assert_eq!(new, patch.observe_waveguide(&n, 1).unwrap()); // reading is not stepping
            let node = n.node_frame(1).unwrap();
            let vm = 0.5 * (old.velocity_m_s + new.velocity_m_s);
            let xm = 0.5 * (old.displacement_m + new.displacement_m);
            let inertia = patch.wall.surface_density * (new.velocity_m_s - old.velocity_m_s) / dt;
            let drag = patch.wall.resistance * vm;
            let spring = patch.wall.stiffness_per_area * xm;
            let scale =
                (inertia.abs() + drag.abs() + spring.abs() + node.pressure_pa.abs()).max(1e-30);
            assert!((inertia + drag + spring - node.pressure_pa).abs() < 2e-10 * scale);
            let displacement_scale =
                (new.displacement_m.abs() + old.displacement_m.abs() + (dt * vm).abs()).max(1e-30);
            assert!(
                (new.displacement_m - old.displacement_m - dt * vm).abs()
                    < 2e-11 * displacement_scale
            );
            close(node.load_flow_m3_s, patch.area_m2 * vm);
            close(node.stored_energy_j, new.stored_energy_j);
            close(
                node.absorbed_energy_j,
                patch.area_m2 * patch.wall.resistance * vm * vm * dt,
            );
            assert!(f.balance_residual_j().abs() < 1e-18);
            peak = peak.max(new.displacement_m.abs());
            total_loss += node.absorbed_energy_j;
            trace.push(n.node_frame(2).unwrap().net_flow_into_node_m3_s);
            old = new;
        }
        assert!(peak > 0.0 && peak < 0.001 * walls[0].radius_m());
        assert!(total_loss > 0.0);
        signals.push(trace);
    }
    assert_ne!(signals[0], signals[1]); // only the material E changed
}

#[path = "material_duct_wall/hereditary.rs"]
mod hereditary;

#[path = "material_duct_wall/spatial.rs"]
mod spatial;

#[path = "material_duct_wall/reeds.rs"]
mod reeds;

/// G1/G3: changing only source E changes both the real pressure trajectory and
/// the wall motion in the original time-domain DG solver. Its closed uniform
/// mode has an independent radial oscillator solution, not a fitted sound.
#[test]
fn source_material_wall_drives_native_dg_pressure_and_storage() {
    use fs_couple::bernoulli_aperture::wall::material::LocalMaterialWall;
    use fs_duct::dg::{DgBoundary, DgDuct, DgMesh};
    let radius = 0.006;
    let length = 0.3;
    let geometry = fs_duct::Duct {
        segments: vec![fs_duct::Segment::Cylinder { radius, length }],
    };
    let gas = gas();
    let mut histories = Vec::new();
    for young in [1e6, 4e6] {
        let source = state(&properties(1000.0, young, 0.0), None, None);
        let wall = LocalMaterialWall::new(
            source.clone(),
            0.0002,
            AxialConstraint::StressFree,
            WallMaterialLoss::KelvinVoigt,
            BAND,
        )
        .unwrap();
        let dg = DgDuct::new(
            &geometry,
            &gas,
            DgMesh {
                elements: 2,
                order: 2,
            },
            DgBoundary::Closed,
            DgBoundary::Closed,
        )
        .unwrap();
        let bound = wall.bind_dg(dg).unwrap();
        assert_eq!(bound.wall().material().identity(), source.identity());
        let dg = bound.duct();
        let pin = dg.compliant_wall().unwrap();
        close(pin.stiffness_per_area, young * 0.0002 / radius.powi(2));
        close(pin.surface_density, 0.2);
        let beta = 2.0 * gas.density * gas.sound_speed.powi(2) / radius;
        let omega = ((pin.stiffness_per_area + beta) / pin.surface_density).sqrt();
        let dt = 1e-6;
        let propagator = dg.propagator(dt, 1e-13).unwrap();
        let mut current = dg.project_initial_state(|_| 1.0, |_| 0.0);
        let initial = current.stored_energy_j();
        let mut pressure = Vec::new();
        for step in 1..=100 {
            propagator.step(&mut current, &[]).unwrap();
            propagator.sync_nodal(&mut current).unwrap();
            let w = (1.0 - (omega * step as f64 * dt).cos()) / (pin.stiffness_per_area + beta);
            assert!((current.pressure()[0] - (1.0 - beta * w)).abs() < 1e-10);
            close(current.stored_energy_j(), initial);
            let wall = dg.compliant_wall_state(&current).unwrap();
            assert!((wall.displacement_m[0] - w).abs() < 1e-17);
            pressure.push(current.pressure()[0]);
        }
        assert!(dg.compliant_wall_state(&current).unwrap().stored_energy_j > 0.0);
        histories.push(pressure);
    }
    assert_ne!(histories[0], histories[1]);
}
