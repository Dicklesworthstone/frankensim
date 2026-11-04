//! V-13a/V-13b battery (bead wf-root-guzez.5.13.1, E4.6a-i): arithmetic +
//! work closure per ITEM and total, the permutation-injection proof that
//! the per-item oracles are live, RSS uncertainty vs hand computation,
//! admission caps at cap AND cap+1, the explicit interference line, the
//! flat-plate fallback as a distinct mode, and the V-13a power-balance
//! band at the Dec-17 state. Repro:
//! cargo test -p fs-flyer --test dragledger_battery

use fs_flyer::dragledger::{
    ComponentDragGeometry, ComponentDragShape, DragComponent, DragLedger, MAX_COMPONENTS,
    flat_plate_aggregate, wright_ledger_v1,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v13\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294; // Dec-17 cold sea-level (air-state-v1)
const V: f64 = 13.86; // airspeed_dec17_mps

fn spatial_ledger() -> DragLedger {
    DragLedger {
        components: ["first", "second"]
            .map(|component_id| DragComponent {
                component_id,
                area_m2: 2.0,
                cd: 1.0,
                cd_source: "declared numerical pressure-face fixture",
                uncertainty_frac: 0.2,
            })
            .to_vec(),
        unresolved_interference_frac: 0.1,
    }
}

#[test]
fn g0_directional_faces_use_normal_speed_and_close_spatial_work() {
    let mut ledger = spatial_ledger();
    ledger.components.truncate(1);
    ledger.unresolved_interference_frac = 0.25;
    let mut geometry = [ComponentDragGeometry {
        component_id: "first",
        point_body_m: [1.0, -2.0, 0.5],
        shape: ComponentDragShape::OrthogonalFaces {
            normals_body: [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            areas_m2: [2.0, 0.5, 0.0],
        },
    }];
    let flow = [10.0, -3.0, 4.0];
    let load = ledger
        .wrench_with_geometry(2.0, &geometry, [0.0; 3], &|_| Ok(flow))
        .unwrap();
    // Independent two-sided pressure faces: Fy=+2*3², Fz=-0.5*4²;
    // the explicit 25% interference allowance scales force and work together.
    assert_eq!(load.force_n, [0.0, 22.5, -10.0]);
    assert_eq!(load.moment_nm, [8.75, 10.0, 22.5]);
    assert_eq!(load.dissipated_power_w, 107.5);
    assert!((load.drag_n - 107.5 / 125.0_f64.sqrt()).abs() < 1e-12);
    let tangent = ledger
        .wrench_with_geometry(2.0, &geometry, [0.0; 3], &|_| Ok([100.0, -3.0, 4.0]))
        .unwrap();
    assert_eq!(
        tangent.force_n, load.force_n,
        "tangential speed cannot invent normal pressure"
    );
    assert_eq!(tangent.dissipated_power_w, load.dissipated_power_w);
    for velocity in [[0.0; 3], [100.0, 0.0, 0.0]] {
        assert_eq!(
            ledger
                .wrench_with_geometry(2.0, &geometry, [0.0; 3], &|_| Ok(velocity))
                .unwrap(),
            Default::default()
        );
    }
    let reverse = ledger
        .wrench_with_geometry(2.0, &geometry, [0.0; 3], &|_| Ok(flow.map(|v| -v)))
        .unwrap();
    assert_eq!(reverse.force_n, load.force_n.map(|v| -v));
    assert_eq!(reverse.moment_nm, load.moment_nm.map(|v| -v));
    assert_eq!(reverse.dissipated_power_w, load.dissipated_power_w);
    let dot = |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let velocity = [1.0, 2.0, 3.0];
    let omega = [0.2, -0.1, 0.3];
    let point_velocity = [1.55, 2.2, 2.7]; // V + omega cross [1,-2,0.5]
    let wind = std::array::from_fn(|i| point_velocity[i] - flow[i]);
    let work = dot(load.force_n, velocity) + dot(load.moment_nm, omega) - dot(load.force_n, wind);
    assert!((work + load.dissipated_power_w).abs() < 1e-12);
    // A proper cyclic rotation of the entire physical fixture rotates its
    // force and moment, while relative loss remains a frame scalar.
    let rotate = |v: [f64; 3]| [v[1], v[2], v[0]];
    geometry[0].point_body_m = rotate(geometry[0].point_body_m);
    if let ComponentDragShape::OrthogonalFaces { normals_body, .. } = &mut geometry[0].shape {
        *normals_body = normals_body.map(rotate);
    }
    let rotated = ledger
        .wrench_with_geometry(2.0, &geometry, [0.0; 3], &|_| Ok(rotate(flow)))
        .unwrap();
    assert_eq!(rotated.force_n, rotate(load.force_n));
    assert_eq!(rotated.moment_nm, rotate(load.moment_nm));
    assert_eq!(rotated.dissipated_power_w, load.dissipated_power_w);
}

#[test]
fn g0_separated_components_resist_rotation_without_inventing_net_translation() {
    let ledger = spatial_ledger();
    let geometry = [
        ComponentDragGeometry {
            component_id: "first",
            point_body_m: [1.5, 0.0, 0.0],
            shape: ComponentDragShape::Isotropic,
        },
        ComponentDragGeometry {
            component_id: "second",
            point_body_m: [-1.5, 0.0, 0.0],
            shape: ComponentDragShape::Isotropic,
        },
    ];
    // Pure yaw omega=2 rad/s in still air. Each component sees opposite
    // crossflow of magnitude 3, even though the CG has zero airspeed.
    let load = ledger
        .wrench_with_geometry(2.0, &geometry, [0.0; 3], &|point| {
            Ok([0.0, 2.0 * point[0], 0.0])
        })
        .unwrap();
    assert_eq!(load.force_n, [0.0; 3]);
    assert_eq!(load.moment_nm[0], 0.0);
    assert_eq!(load.moment_nm[1], 0.0);
    assert!((load.moment_nm[2] + 59.4).abs() < 1e-12);
    assert!((load.dissipated_power_w - 118.8).abs() < 1e-12);
    assert!((2.0 * load.moment_nm[2] + load.dissipated_power_w).abs() < 1e-12);
}

#[test]
fn g0_component_geometry_rejects_omission_rebinding_and_invalid_faces_before_sampling() {
    let ledger = spatial_ledger();
    let valid = [
        ComponentDragGeometry {
            component_id: "first",
            point_body_m: [1.0, 0.0, 0.0],
            shape: ComponentDragShape::Isotropic,
        },
        ComponentDragGeometry {
            component_id: "second",
            point_body_m: [-1.0, 0.0, 0.0],
            shape: ComponentDragShape::Isotropic,
        },
    ];
    let mut invalid = vec![
        valid[..1].to_vec(),
        vec![valid[1], valid[0]],
        vec![valid[0], valid[0]],
    ];
    for shape in [
        ComponentDragShape::OrthogonalFaces {
            normals_body: [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            areas_m2: [1.0; 3],
        },
        ComponentDragShape::OrthogonalFaces {
            normals_body: [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            areas_m2: [1.0; 3],
        },
        ComponentDragShape::OrthogonalFaces {
            normals_body: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            areas_m2: [-1.0, 1.0, 1.0],
        },
        ComponentDragShape::OrthogonalFaces {
            normals_body: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            areas_m2: [0.0; 3],
        },
    ] {
        let mut bad = valid;
        bad[0].shape = shape;
        invalid.push(bad.to_vec());
    }
    let mut bad_point = valid;
    bad_point[0].point_body_m[2] = f64::NAN;
    invalid.push(bad_point.to_vec());
    for geometry in invalid {
        let result = ledger.wrench_with_geometry(2.0, &geometry, [0.0; 3], &|_| {
            panic!("invalid geometry sampled the air")
        });
        assert_eq!(result.unwrap_err().code, "drag-geometry-invalid");
    }
    let refused = ledger.wrench_with_geometry(2.0, &valid, [0.0; 3], &|_| {
        Err(fs_flyer::Refusal {
            code: "component-wind-refused",
            message: "fixture".into(),
            ranked_repairs: vec![],
        })
    });
    assert_eq!(refused.unwrap_err().code, "component-wind-refused");
    assert_eq!(
        ledger
            .wrench_with_geometry(2.0, &valid, [0.0; 3], &|_| Ok([f64::NAN, 0.0, 0.0]))
            .unwrap_err()
            .code,
        "air-state-invalid"
    );
    assert_eq!(
        ledger
            .wrench_with_geometry(2.0, &valid, [0.0; 3], &|_| Ok([1e200, 0.0, 0.0]))
            .unwrap_err()
            .code,
        "drag-wrench-overflow"
    );
}

#[test]
fn g0_resolved_isotropic_layout_retains_uniform_ledger_and_separate_geometry_identity() {
    let ledger = wright_ledger_v1();
    let original_digest = ledger.digest();
    let mut geometry: Vec<_> = ledger
        .components
        .iter()
        .map(|c| ComponentDragGeometry {
            component_id: c.component_id,
            point_body_m: [0.0, 0.0, 0.4],
            shape: ComponentDragShape::Isotropic,
        })
        .collect();
    let velocity = [3.0, 4.0, 12.0];
    let resolved = ledger
        .wrench_with_geometry(RHO, &geometry, [0.0; 3], &|_| Ok(velocity))
        .unwrap();
    let aggregate = ledger.wrench_at(RHO, velocity, [0.0, 0.0, 0.4]).unwrap();
    for (a, b) in resolved
        .force_n
        .into_iter()
        .chain(resolved.moment_nm)
        .chain([resolved.drag_n, resolved.dissipated_power_w])
        .zip(
            aggregate
                .force_n
                .into_iter()
                .chain(aggregate.moment_nm)
                .chain([aggregate.drag_n, aggregate.dissipated_power_w]),
        )
    {
        assert!((a - b).abs() < 1e-11);
    }
    let spatial = ledger.spatial_digest(&geometry);
    geometry[0].point_body_m[0] = 0.2;
    assert_ne!(spatial, ledger.spatial_digest(&geometry));
    let moved = ledger.spatial_digest(&geometry);
    geometry[0].shape = ComponentDragShape::OrthogonalFaces {
        normals_body: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        areas_m2: [0.3, 0.1, 0.2],
    };
    assert_ne!(moved, ledger.spatial_digest(&geometry));
    assert_eq!(ledger.digest(), original_digest);
}

#[test]
fn g0_local_drag_wrench_obeys_spatial_work_and_reverses_with_relative_flow() {
    let ledger = wright_ledger_v1();
    let flow = [3.0, 4.0, 12.0];
    let arm = [0.2, -0.4, 0.6];
    let wrench = ledger.wrench_at(RHO, flow, arm).unwrap();
    let report = ledger.evaluate(RHO, 13.0).unwrap();
    assert_eq!(wrench.drag_n, report.total_parasite_n);
    assert_eq!(wrench.dissipated_power_w, report.power_w);
    let dot = |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    assert!((dot(wrench.force_n, flow) + report.power_w).abs() < 1e-9);
    // V_point = V_origin + omega x arm, with the material wind subtracted.
    let velocity = [10.0, -2.0, 3.0];
    let omega = [0.7, -0.3, 0.4];
    let point = [
        velocity[0] + omega[1] * arm[2] - omega[2] * arm[1],
        velocity[1] + omega[2] * arm[0] - omega[0] * arm[2],
        velocity[2] + omega[0] * arm[1] - omega[1] * arm[0],
    ];
    let wind = std::array::from_fn(|i| point[i] - flow[i]);
    let body_power = dot(wrench.force_n, velocity) + dot(wrench.moment_nm, omega);
    let wind_power = dot(wrench.force_n, wind);
    assert!((body_power - wind_power + wrench.dissipated_power_w).abs() < 1e-9);
    let reversed = ledger.wrench_at(RHO, flow.map(|v| -v), arm).unwrap();
    assert_eq!(reversed.force_n, wrench.force_n.map(|v| -v));
    assert_eq!(reversed.moment_nm, wrench.moment_nm.map(|v| -v));
    assert_eq!(reversed.dissipated_power_w, wrench.dissipated_power_w);
    assert_eq!(
        ledger.wrench_at(RHO, [0.0; 3], arm).unwrap(),
        Default::default()
    );
}

#[test]
fn g0_local_drag_retains_representable_loads_at_extreme_velocity_scales() {
    let ledger = wright_ledger_v1();
    for (rho, speed) in [(1e300, 1e-200), (1e-300, 1e200)] {
        let wrench = ledger.wrench_at(rho, [speed, 0.0, 0.0], [0.0; 3]).unwrap();
        let report = ledger.evaluate(rho, speed).unwrap();
        assert!(wrench.drag_n > 0.0);
        assert_eq!(wrench.drag_n, report.total_parasite_n);
        assert_eq!(wrench.force_n[0], -report.total_parasite_n);
        assert_eq!(wrench.dissipated_power_w, report.power_w);
    }
}

#[test]
fn g0_local_drag_refuses_invalid_and_overflowing_loads() {
    let ledger = wright_ledger_v1();
    for (rho, flow, arm) in [
        (0.0, [1.0; 3], [0.0; 3]),
        (RHO, [f64::NAN, 0.0, 0.0], [0.0; 3]),
        (RHO, [1.0; 3], [0.0, f64::INFINITY, 0.0]),
    ] {
        assert_eq!(
            ledger.wrench_at(rho, flow, arm).unwrap_err().code,
            "air-state-invalid"
        );
    }
    assert_eq!(
        ledger
            .wrench_at(RHO, [1e200, 0.0, 0.0], [0.0; 3])
            .unwrap_err()
            .code,
        "drag-wrench-overflow"
    );
}

#[test]
fn per_item_arithmetic_and_work_closure() {
    let rep = wright_ledger_v1().evaluate(RHO, V).unwrap();
    let q = 0.5 * RHO * V * V;
    let ledger = wright_ledger_v1();
    // PER-ITEM oracles (never totals-only): each line is q*A*cd and
    // P = D*V, checked against an independent hand computation.
    for (item, comp) in rep.items.iter().zip(&ledger.components) {
        assert_eq!(item.component_id, comp.component_id);
        let d_ref = q * comp.area_m2 * comp.cd;
        assert!(
            (item.drag_n - d_ref).abs() < 1e-12 * d_ref,
            "{}: {} vs {d_ref}",
            comp.component_id,
            item.drag_n
        );
        assert!(
            (item.power_w - item.drag_n * V).abs() < 1e-9,
            "work closure per item"
        );
        assert!((item.sigma_n - d_ref * comp.uncertainty_frac).abs() < 1e-9);
    }
    // Total closure: components + the EXPLICIT interference line.
    let sum: f64 = rep.items.iter().map(|i| i.drag_n).sum();
    assert!((rep.component_sum_n - sum).abs() < 1e-9);
    assert!(
        (rep.unresolved_interference_drag_n - sum * 0.10).abs() < 1e-9,
        "interference must be its own line at the declared fraction"
    );
    assert!((rep.total_parasite_n - (sum + rep.unresolved_interference_drag_n)).abs() < 1e-9);
    assert!((rep.power_w - rep.total_parasite_n * V).abs() < 1e-6);
    // RSS by hand.
    let mut var: f64 = rep.items.iter().map(|i| i.sigma_n * i.sigma_n).sum();
    var += rep.unresolved_interference_drag_n * rep.unresolved_interference_drag_n;
    assert!((rep.sigma_total_n - var.sqrt()).abs() < 1e-9);
    jlog(
        "closure",
        &format!(
            "\"total_n\":{},\"interference_n\":{},\"sigma_n\":{},\"power_w\":{}",
            rep.total_parasite_n,
            rep.unresolved_interference_drag_n,
            rep.sigma_total_n,
            rep.power_w
        ),
    );
}

#[test]
fn permutation_injection_proves_item_oracles_live() {
    // Sum tests are blind to permutation (memory doctrine): swap two
    // components' AREAS and verify the per-item oracle CATCHES it while
    // the total stays identical — proof the per-item checks are the load
    // bearers, executed, not asserted.
    let mut swapped = wright_ledger_v1();
    let (a0, a1) = (swapped.components[0].area_m2, swapped.components[1].area_m2);
    let (c0, c1) = (swapped.components[0].cd, swapped.components[1].cd);
    // Swap the (area*cd) PRODUCTS between lines 0 and 1.
    swapped.components[0].area_m2 = a1;
    swapped.components[0].cd = c1;
    swapped.components[1].area_m2 = a0;
    swapped.components[1].cd = c0;
    let clean = wright_ledger_v1().evaluate(RHO, V).unwrap();
    let bad = swapped.evaluate(RHO, V).unwrap();
    assert!(
        (clean.component_sum_n - bad.component_sum_n).abs() < 1e-9,
        "the injected swap must be INVISIBLE to the total"
    );
    let mut caught = 0;
    for (ci, bi) in clean.items.iter().zip(&bad.items) {
        if (ci.drag_n - bi.drag_n).abs() > 1e-9 {
            caught += 1;
        }
    }
    assert!(
        caught >= 2,
        "per-item oracles must catch the swap the total missed"
    );
    jlog("permutation", &format!("\"caught_items\":{caught}"));
}

#[test]
fn admission_caps_at_cap_and_cap_plus_one() {
    let mk = |n: usize| -> DragLedger {
        // Distinct ids from a static table (component_id is &'static str).
        const IDS: [&str; 33] = [
            "c00", "c01", "c02", "c03", "c04", "c05", "c06", "c07", "c08", "c09", "c10", "c11",
            "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22", "c23",
            "c24", "c25", "c26", "c27", "c28", "c29", "c30", "c31", "c32",
        ];
        DragLedger {
            components: (0..n)
                .map(|i| DragComponent {
                    component_id: IDS[i],
                    area_m2: 0.1,
                    cd: 1.0,
                    cd_source: "test",
                    uncertainty_frac: 0.2,
                })
                .collect(),
            unresolved_interference_frac: 0.1,
        }
    };
    assert!(mk(MAX_COMPONENTS).admit().is_ok(), "cap admits");
    assert_eq!(
        mk(MAX_COMPONENTS + 1).admit().unwrap_err().code,
        "ledger-component-count-invalid",
        "cap+1 refuses"
    );
    assert_eq!(
        mk(0).admit().unwrap_err().code,
        "ledger-component-count-invalid",
        "empty refuses"
    );
    // Component field refusals.
    let mut bad = wright_ledger_v1();
    bad.components[0].area_m2 = -0.1;
    assert_eq!(
        bad.evaluate(RHO, V).unwrap_err().code,
        "ledger-component-invalid"
    );
    let mut bad = wright_ledger_v1();
    bad.components[0].uncertainty_frac = 1.0 + f64::EPSILON;
    assert_eq!(bad.admit().unwrap_err().code, "ledger-component-invalid");
    let mut bad = wright_ledger_v1();
    bad.components[0].uncertainty_frac = 1.0;
    assert!(bad.admit().is_ok(), "uncertainty exactly 1 is the cap");
    // Duplicate id.
    let mut dup = wright_ledger_v1();
    dup.components[1].component_id = dup.components[0].component_id;
    assert_eq!(dup.admit().unwrap_err().code, "ledger-component-duplicate");
    // Interference window at cap and past it.
    let mut edge = wright_ledger_v1();
    edge.unresolved_interference_frac = 0.5;
    assert!(edge.admit().is_ok());
    edge.unresolved_interference_frac = 0.5 + f64::EPSILON;
    assert_eq!(
        edge.admit().unwrap_err().code,
        "ledger-interference-invalid"
    );
    // Air state.
    assert_eq!(
        wright_ledger_v1().evaluate(RHO, 0.0).unwrap_err().code,
        "air-state-invalid"
    );
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn flat_plate_fallback_is_a_distinct_identified_mode() {
    // The fallback returns a DIFFERENT type with the equivalent area
    // identified — and for the same equivalent area it reproduces the
    // ledger total (consistency), while carrying no allocation claims.
    let rep = wright_ledger_v1().evaluate(RHO, V).unwrap();
    let q = 0.5 * RHO * V * V;
    let f_equiv = rep.total_parasite_n / q;
    let fb = flat_plate_aggregate(f_equiv, RHO, V).unwrap();
    assert!((fb.drag_n - rep.total_parasite_n).abs() < 1e-9);
    assert!((fb.flat_plate_area_m2 - f_equiv).abs() < 1e-12);
    assert_eq!(
        flat_plate_aggregate(0.0, RHO, V).unwrap_err().code,
        "flat-plate-area-invalid"
    );
    jlog("fallback", &format!("\"f_equiv_m2\":{f_equiv}"));
}

#[test]
fn v13a_power_balance_band_at_dec17() {
    // V-13a: the INDEPENDENT full-aircraft evidence is the power balance
    // of a barely-sustained flight: engine_power_w = 8948 (verified),
    // prop efficiency ~0.66 (Wright bench class), so thrust power
    // available ~5.9 kW at 13.86 m/s. Total drag power required =
    // (parasite ledger + wing profile + induced) * V and must sit inside
    // [4.2, 6.9] kW — a factor-honest band, not fake precision.
    let rep = wright_ledger_v1().evaluate(RHO, V).unwrap();
    let q = 0.5 * RHO * V * V;
    // Wing-owned drag terms, computed from verified geometry (both
    // planes 47.38 m², span 12.29 m -> AR 6.38 biplane) with declared
    // Estimated constants: profile cd 0.015, Oswald-class e 0.7.
    let (s, b) = (47.38, 12.29);
    let ar = b * b / (s / 2.0); // per-plane AR on shared span convention
    let weight_n = 340.2 * 9.80665;
    let cl = weight_n / (q * s);
    let cdi = cl * cl / (core::f64::consts::PI * ar * 0.7);
    let d_wing = q * s * (0.015 + cdi);
    let total_drag = rep.total_parasite_n + d_wing;
    let p_required = total_drag * V;
    assert!(
        (4200.0..=6900.0).contains(&p_required),
        "V-13a: required power {p_required} W outside the power-balance band"
    );
    // Anti-vacuity twin: a ledger 3x heavier must LEAVE the band —
    // proving the band can fail.
    let mut bloated = wright_ledger_v1();
    for c in &mut bloated.components {
        c.area_m2 *= 3.0;
    }
    let p_bloated = bloated.evaluate(RHO, V).unwrap().total_parasite_n * V + d_wing * V;
    assert!(
        !(4200.0..=6900.0).contains(&p_bloated),
        "the band must be falsifiable: bloated ledger scored {p_bloated} W inside it"
    );
    jlog(
        "v13a",
        &format!(
            "\"parasite_n\":{},\"wing_n\":{d_wing},\"p_required_w\":{p_required},\"band\":[4200,6900],\"cl\":{cl}",
            rep.total_parasite_n
        ),
    );
}

#[test]
fn determinism_and_golden_digest() {
    let a = wright_ledger_v1().evaluate(RHO, V).unwrap();
    let b = wright_ledger_v1().evaluate(RHO, V).unwrap();
    assert_eq!(a, b, "bitwise repeat");
    let mut payload = Vec::new();
    for i in &a.items {
        payload.extend_from_slice(&i.drag_n.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(&a.total_parasite_n.to_bits().to_le_bytes());
    payload.extend_from_slice(&a.sigma_total_n.to_bits().to_le_bytes());
    payload.extend_from_slice(a.ledger_digest.as_bytes());
    let digest = fs_blake3::hash_domain("org.frankensim.fs-flyer.v13-golden.v1", &payload).to_hex();
    // Digest sensitivity: a changed component moves the ledger digest.
    let mut other = wright_ledger_v1();
    other.components[0].area_m2 *= 1.5;
    assert_ne!(other.digest(), wright_ledger_v1().digest());
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "c3cbf9b93d49fbcf9ea70955c6734da4008a8ca97647cf83f5509b17b05e3f1d",
        "V-13 ledger golden moved — determinism regression or an \
         intentional ledger change requiring the golden-bump protocol"
    );
}
