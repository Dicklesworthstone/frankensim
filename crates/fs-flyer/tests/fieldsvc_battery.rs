//! E7.1-i battery (bead wf-root-guzez.8.1.1): the §5.5 field service.
//! Exact-vs-FD derivative duals VERIFIED per point per entry on the
//! analytic atmosphere (never a totals-only norm); tensor-invariant
//! oracles (Q recomputed, λ₂ char-poly residual); component
//! additivity bitwise; omitted-component honesty with the
//! forbidden-claims falsifier EXECUTED both ways; singularity-core +
//! normalized-divergence masking; ground-image on-plane cancellation
//! (V-06a echo); caps at cap AND cap+1; determinism golden.
//! Repro: cargo test -p fs-flyer --test fieldsvc_battery

use fs_atmo::{Atmosphere, DEC17_AIR, FlatSiteLogLaw, TurbulenceField};
use fs_flyer::fieldsvc::{
    BoundSystem, C_BOUND_CIRCULATION, C_GROUND_IMAGES, C_GUST_EVENT, C_MEAN_ATMO, C_TURB_ATMO,
    FieldSourceStateV1, GridSpec, MAX_POINTS, claim_total_flow, sample_field,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-fieldsvc\",\"case\":\"{case}\",{payload}}}");
}

fn atmosphere(n_modes: usize) -> Atmosphere {
    Atmosphere {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 5.0e-3,
            displacement_height_m: 0.02,
            reference_height_m: 10.0,
            reference_speed_mps: 8.0,
        },
        turbulence: TurbulenceField::build(1903, n_modes, 0.9, 20.0, 8.0).unwrap(),
        air: DEC17_AIR,
    }
}

fn ambient_state() -> FieldSourceStateV1 {
    FieldSourceStateV1 {
        tick: 240,
        source_state_digest: "battery-source-digest".into(),
        atmosphere: atmosphere(16),
        bound: None,
        images_active: false,
    }
}

fn lifting_state() -> FieldSourceStateV1 {
    FieldSourceStateV1 {
        bound: Some(BoundSystem {
            gamma_m2ps: 8.0,
            tip_left_m: [0.0, -6.0, 2.5],
            tip_right_m: [0.0, 6.0, 2.5],
            trail_m: 40.0,
            core_m: 0.05,
        }),
        images_active: true,
        ..ambient_state()
    }
}

fn grid_small() -> GridSpec {
    GridSpec {
        origin_m: [-2.0, -3.0, 0.8],
        dx_m: 1.1,
        nx: 4,
        ny: 3,
        nz: 3,
    }
}

#[test]
fn exact_vs_fd_dual_verified_per_point() {
    let state = ambient_state();
    let s = sample_field(&state, &grid_small(), C_MEAN_ATMO | C_TURB_ATMO).unwrap();
    let g = grid_small();
    let h = 1.0e-4;
    let mut worst_grad = 0.0f64;
    let mut worst_div = 0.0f64;
    for i in 0..g.n_points() {
        // The analytic construction is solenoidal (curl field + pure
        // shear): the analytic divergence must vanish to roundoff.
        assert!(
            s.div_analytic[i].abs() < 1e-10,
            "point {i}: analytic div {}",
            s.div_analytic[i]
        );
        assert!(
            (s.div_finite_difference[i] - s.div_analytic[i]).abs() < 1e-5,
            "point {i}: FD dual {} vs {}",
            s.div_finite_difference[i],
            s.div_analytic[i]
        );
        // Per-ENTRY gradient dual: independent central differences of
        // the sampled velocity against the analytic grad_u.
        let p = g.point(i);
        for j in 0..3 {
            let mut pp = p;
            let mut pm = p;
            pp[j] += h;
            pm[j] -= h;
            let up = probe_u(&state, pp);
            let um = probe_u(&state, pm);
            for c in 0..3 {
                let fd = (up[c] - um[c]) / (2.0 * h);
                let d = (fd - s.grad_u[i][c][j]).abs();
                worst_grad = worst_grad.max(d);
                assert!(
                    d < 1e-5,
                    "point {i} d u{c}/dx{j}: fd {fd} vs {}",
                    s.grad_u[i][c][j]
                );
            }
        }
        // Vorticity is the antisymmetric part of the SAME gradient.
        let gu = &s.grad_u[i];
        assert_eq!(s.omega[i][0].to_bits(), (gu[2][1] - gu[1][2]).to_bits());
        worst_div = worst_div.max((s.div_finite_difference[i] - s.div_analytic[i]).abs());
    }
    jlog(
        "dual",
        &format!("\"worst_grad_dual\":{worst_grad},\"worst_div_dual\":{worst_div}"),
    );
}

/// Independent velocity probe for the FD cross-check (mean + turb).
fn probe_u(state: &FieldSourceStateV1, p: [f64; 3]) -> [f64; 3] {
    let ts = state
        .atmosphere
        .turbulence
        .sample(p[0], p[1], p[2], state.tick);
    [
        ts.u[0] + state.atmosphere.mean.speed(p[2]),
        ts.u[1],
        ts.u[2],
    ]
}

#[test]
fn tensor_invariants_recomputed_per_point() {
    let s = sample_field(&ambient_state(), &grid_small(), C_MEAN_ATMO | C_TURB_ATMO).unwrap();
    for i in 0..s.u.len() {
        let g = &s.grad_u[i];
        let mut s2 = 0.0;
        let mut w2 = 0.0;
        let mut a = [[0.0f64; 3]; 3];
        let mut sm = [[0.0f64; 3]; 3];
        let mut wm = [[0.0f64; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                sm[r][c] = 0.5 * (g[r][c] + g[c][r]);
                wm[r][c] = 0.5 * (g[r][c] - g[c][r]);
                s2 += sm[r][c] * sm[r][c];
                w2 += wm[r][c] * wm[r][c];
            }
        }
        for r in 0..3 {
            for c in 0..3 {
                for k in 0..3 {
                    a[r][c] += sm[r][k] * sm[k][c] + wm[r][k] * wm[k][c];
                }
            }
        }
        assert_eq!(
            s.q_criterion[i].to_bits(),
            (0.5 * (w2 - s2)).to_bits(),
            "Q at {i}"
        );
        // λ₂ satisfies the characteristic polynomial of A.
        let l = s.lambda2[i];
        let det3 = |m: &[[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let mut al = a;
        for d in 0..3 {
            al[d][d] -= l;
        }
        let scale = a.iter().flatten().map(|v| v.abs()).fold(1e-30, f64::max);
        assert!(
            det3(&al).abs() / (scale * scale * scale) < 1e-9,
            "lambda2 char-poly residual at {i}"
        );
        assert!(s.strain_magnitude[i] >= 0.0);
    }
    jlog("invariants", &format!("\"points\":{}", s.u.len()));
}

#[test]
fn component_additivity_and_claim_honesty() {
    let state = lifting_state();
    let g = grid_small();
    let both = sample_field(&state, &g, C_MEAN_ATMO | C_TURB_ATMO).unwrap();
    let mean = sample_field(&state, &g, C_MEAN_ATMO).unwrap();
    let turb = sample_field(&state, &g, C_TURB_ATMO).unwrap();
    for i in 0..g.n_points() {
        for c in 0..3 {
            assert_eq!(
                both.u[i][c].to_bits(),
                (mean.u[i][c] + turb.u[i][c]).to_bits(),
                "additivity at {i}.{c}"
            );
        }
    }
    // Honesty: the state SUPPORTS bound circulation; a sum without it
    // names the omission and the total-flow claim REFUSES.
    assert!(both.meta.omitted_components.contains(&"bound-circulation"));
    let err = claim_total_flow(&both).unwrap_err();
    assert_eq!(err.code, "forbidden-claim-total-flow");
    assert!(err.message.contains("bound-circulation"));
    // Including every SUPPORTED force-coupled component allows it —
    // unsupported branches (physical wake, prop) do not forbid.
    let full = sample_field(
        &state,
        &g,
        C_MEAN_ATMO | C_TURB_ATMO | C_BOUND_CIRCULATION | C_GROUND_IMAGES,
    )
    .unwrap();
    claim_total_flow(&full).unwrap();
    assert!(full.meta.omitted_components.contains(&"physical-wake"));
    assert_eq!(full.provenance.len(), 4);
    jlog(
        "honesty",
        &format!(
            "\"refusal\":\"{}\",\"omitted_full\":{}",
            err.code,
            full.meta.omitted_components.len()
        ),
    );
}

#[test]
fn singular_core_masking_and_image_cancellation() {
    let state = lifting_state();
    // A vertical line of points crossing the span at mid-wing.
    let g = GridSpec {
        origin_m: [0.0, 0.0, 2.30],
        dx_m: 0.05,
        nx: 1,
        ny: 1,
        nz: 9,
    };
    let s = sample_field(&state, &g, C_MEAN_ATMO | C_BOUND_CIRCULATION).unwrap();
    let mut cores = 0;
    for i in 0..9 {
        let z = 2.30 + 0.05 * i as f64;
        let near = (z - 2.5).abs() < 4.0 * 0.05;
        assert_eq!(s.singularity_core_mask[i], near, "core mask at z {z}");
        if near {
            cores += 1;
            assert!(
                s.normalized_divergence(i).is_none(),
                "normalized div must be masked in core"
            );
        }
    }
    assert!(cores >= 3, "the line must actually cross the core");
    assert_eq!(s.meta.core_radius_m, 0.05);
    // Ground images: on the certified plane the normal velocity of
    // bound + images cancels (V-06a doctrine echoed by this service).
    let gp = GridSpec {
        origin_m: [-3.0, -4.0, 0.0],
        dx_m: 2.5,
        nx: 3,
        ny: 4,
        nz: 1,
    };
    let si = sample_field(&state, &gp, C_BOUND_CIRCULATION | C_GROUND_IMAGES).unwrap();
    for i in 0..12 {
        let speed = (si.u[i][0].powi(2) + si.u[i][1].powi(2) + si.u[i][2].powi(2)).sqrt();
        assert!(
            si.u[i][2].abs() <= 1e-9 * speed.max(1e-12),
            "on-plane normal residual at {i}: {} of {speed}",
            si.u[i][2]
        );
    }
    jlog("core-and-images", &format!("\"core_points\":{cores}"));
}

#[test]
fn refusals_at_cap_and_beyond() {
    let state = ambient_state();
    // AT the point cap admits (mean-only keeps it cheap)…
    let at_cap = GridSpec {
        origin_m: [0.0, 0.0, 1.0],
        dx_m: 0.001,
        nx: MAX_POINTS,
        ny: 1,
        nz: 1,
    };
    assert!(sample_field(&state, &at_cap, C_MEAN_ATMO).is_ok());
    // …one more point refuses.
    let over = GridSpec {
        nx: MAX_POINTS + 1,
        ..at_cap
    };
    assert_eq!(
        sample_field(&state, &over, C_MEAN_ATMO).unwrap_err().code,
        "grid-spec-invalid"
    );
    for bad in [
        GridSpec {
            dx_m: 0.0,
            ..grid_small()
        },
        GridSpec {
            origin_m: [f64::NAN, 0.0, 1.0],
            ..grid_small()
        },
        GridSpec {
            nx: 0,
            ..grid_small()
        },
    ] {
        assert_eq!(
            sample_field(&state, &bad, C_MEAN_ATMO).unwrap_err().code,
            "grid-spec-invalid"
        );
    }
    assert_eq!(
        sample_field(&state, &grid_small(), 0).unwrap_err().code,
        "component-mask-empty"
    );
    assert_eq!(
        sample_field(&state, &grid_small(), 1 << 8)
            .unwrap_err()
            .code,
        "component-mask-unknown"
    );
    // Requesting a component the state cannot produce refuses TYPED —
    // never a silent zero-field.
    let err = sample_field(&state, &grid_small(), C_MEAN_ATMO | C_GUST_EVENT).unwrap_err();
    assert_eq!(err.code, "component-unsupported");
    assert!(err.message.contains("gust-event"));
    // Ambient state without a bound system refuses the bound bit too.
    assert_eq!(
        sample_field(&state, &grid_small(), C_BOUND_CIRCULATION)
            .unwrap_err()
            .code,
        "component-unsupported"
    );
    jlog("refusals", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn solid_exclusion_and_determinism_golden() {
    let state = lifting_state();
    // Straddle the plane: below-plane points are solid-excluded and
    // invalid, never silently zero-valid.
    let g = GridSpec {
        origin_m: [1.0, 1.0, -0.6],
        dx_m: 0.6,
        nx: 2,
        ny: 1,
        nz: 3,
    };
    let s = sample_field(&state, &g, C_MEAN_ATMO | C_BOUND_CIRCULATION).unwrap();
    for i in 0..6 {
        let below = g.point(i)[2] < 0.0;
        assert_eq!(s.solid_exclusion_mask[i], below, "solid at {i}");
        if below {
            assert!(!s.validity_mask[i]);
        }
    }
    // Determinism golden over the full payload.
    let run = || {
        sample_field(
            &state,
            &grid_small(),
            C_MEAN_ATMO | C_TURB_ATMO | C_BOUND_CIRCULATION | C_GROUND_IMAGES,
        )
        .unwrap()
        .digest()
    };
    let a = run();
    assert_eq!(a, run(), "bit-identical twice");
    // Meta identity: id is stable and moves with the source digest.
    let id1 = state.field_source_snapshot_id();
    let mut other = state.clone();
    other.source_state_digest = "another-digest".into();
    assert_ne!(id1, other.field_source_snapshot_id());
    jlog(
        "golden",
        &format!("\"digest\":\"{a}\",\"snapshot_id\":\"{id1}\""),
    );
    assert_eq!(
        a, "e07ab5cb99f09d73349dafbce22510e63177420f6cd7c65433d7b1f70586bdc3",
        "fieldsvc golden moved — determinism regression or an \
         intentional kernel change requiring the golden-bump protocol"
    );
}
