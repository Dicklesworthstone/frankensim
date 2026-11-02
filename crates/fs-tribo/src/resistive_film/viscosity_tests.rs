use super::*;

fn film() -> ResistiveFilm {
    let gap = GapPort {
        reference_m: 0.0002,
        closure: vec![1.0, -1.0],
    };
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 0.0001,
            gap: gap.clone(),
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 0.01,
            length_m: 0.02,
            gap,
        }],
        2,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 0.00001,
            maximum_gap_m: 0.0005,
            maximum_pressure_pa: 10000.0,
        },
    )
    .unwrap()
}

#[test]
fn changed_viscosity_drives_the_same_pressure_force_work_and_tangent() {
    let old = film();
    let next = old.with_viscosity(0.004).unwrap();
    let q = [0.00001, 0.0];
    let v = [0.001, 0.0];
    let mut old_f = [0.0; 2];
    let mut new_f = [0.0; 2];
    let mut old_p = [0.0];
    let mut new_p = [0.0];
    let a = old.evaluate_into(&q, &v, &mut old_f, &mut old_p).unwrap();
    let b = next.evaluate_into(&q, &v, &mut new_f, &mut new_p).unwrap();
    // Independent Poiseuille balance with the actual, already-closed aperture.
    let g = 0.01 * 0.00019_f64.powi(3) / (12.0 * 0.004 * 0.02);
    let p = 0.0001 * 0.001 / g;
    assert!((new_p[0] - p).abs() < 1e-10);
    assert!((new_f[0] - 0.0001 * p).abs() < 1e-14);
    assert_eq!(new_f[1], -new_f[0]);
    assert!((b.dissipated_power_w - 0.4 * a.dissipated_power_w).abs() < 1e-17);
    assert_eq!(a.minimum_gap_m, b.minimum_gap_m);
    let mut da = [0.0; 2];
    let mut db = [0.0; 2];
    old.tangent_into(&q, &v, &[0.000001, 0.0], &[0.0001, 0.0], &mut da)
        .unwrap();
    next.tangent_into(&q, &v, &[0.000001, 0.0], &[0.0001, 0.0], &mut db)
        .unwrap();
    for i in 0..2 {
        assert!((db[i] - 0.4 * da[i]).abs() < 1e-14);
    }
}

#[test]
fn coefficient_recompilation_is_nonaccumulating_and_keeps_the_original_limits() {
    let original = film();
    let mut refreshed = original.clone();
    for _ in 0..64 {
        refreshed = refreshed
            .with_viscosity(0.013)
            .unwrap()
            .with_viscosity(0.01)
            .unwrap();
    }
    assert_eq!(refreshed.viscosity_pa_s(), original.viscosity_pa_s());
    assert_eq!(refreshed.factors, original.factors);
    let mut f = [42.0; 2];
    let mut p = [42.0];
    assert!(
        refreshed
            .evaluate_into(&[0.0002, 0.0], &[0.0; 2], &mut f, &mut p)
            .is_err()
    );
    assert_eq!(f, [42.0; 2]);
    assert_eq!(p, [42.0]);
    assert!(
        refreshed
            .evaluate_into(&[0.0; 2], &[1.0, 0.0], &mut f, &mut p)
            .is_err()
    );
    assert_eq!(f, [42.0; 2]);
    assert_eq!(p, [42.0]);
    for mu in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(original.with_viscosity(mu).is_err());
        assert_eq!(original.viscosity_pa_s(), Some(0.01));
    }
}

// Synthetic coefficients magnify pressure sensitivity; these are exact reduced-
// model oracles, not a calibration or a high-pressure material-data claim.
#[test]
fn barus_single_passage_matches_integrated_fluidity_and_actual_pressure_work() {
    let original = film();
    let alpha = 2e-4;
    let barus = original.with_barus_viscosity(alpha).unwrap();
    for mu in [0.01, 0.004] {
        let current = barus.with_viscosity(mu).unwrap();
        assert_eq!(current.pressure_viscosity_coefficient_per_pa(), alpha);
        assert_eq!(current.viscosity_pa_s(), Some(mu));
        for speed in [-0.001, 0.001] {
            let q = [0.00001, 0.0];
            let mut force = [0.0; 2];
            let mut pressure = [0.0];
            let report = current
                .evaluate_into(&q, &[speed, 0.0], &mut force, &mut pressure)
                .unwrap();
            let g0 = 0.01 * 0.00019_f64.powi(3) / (12.0 * mu * 0.02);
            let volume_rate = 0.0001 * speed;
            // Integrate dp/mu(p) along the channel independently of the
            // pressure solver, including the negative-gauge suction branch.
            let integrated_rate = g0 * (1.0 - (-alpha * pressure[0]).exp()) / alpha;
            assert!((integrated_rate - volume_rate).abs() < 2e-22);
            let exact_pressure = -(1.0 - alpha * volume_rate / g0).ln() / alpha;
            assert!((pressure[0] - exact_pressure).abs() < 2e-11);
            assert!((force[0] - 0.0001 * exact_pressure).abs() < 3e-15);
            assert_eq!(force[1], -force[0]);
            assert!((report.dissipated_power_w - volume_rate * exact_pressure).abs() < 3e-18);
            assert!((report.dissipated_power_w - speed * force[0]).abs() < 3e-18);
            assert!(report.dissipated_power_w > 0.0);
        }
    }
}

fn two_cell_barus() -> ResistiveFilm {
    let gap = |reference_m, closure| GapPort {
        reference_m,
        closure,
    };
    ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: 0.0001,
                gap: gap(0.0002, vec![1.0, 0.0]),
            },
            FilmCell {
                area_m2: 0.0002,
                gap: gap(0.00024, vec![0.0, 1.0]),
            },
        ],
        vec![
            FilmChannel {
                from: 0,
                to: Some(1),
                width_m: 0.01,
                length_m: 0.02,
                gap: gap(0.0002, vec![0.4, 0.6]),
            },
            FilmChannel {
                from: 1,
                to: None,
                width_m: 0.01,
                length_m: 0.02,
                gap: gap(0.00022, vec![0.2, 0.8]),
            },
            FilmChannel {
                from: 0,
                to: None,
                width_m: 0.01,
                length_m: 0.02,
                gap: gap(0.00018, vec![0.7, 0.3]),
            },
        ],
        2,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 0.00001,
            maximum_gap_m: 0.0005,
            maximum_pressure_pa: 10000.0,
        },
    )
    .unwrap()
    .with_barus_viscosity(2e-4)
    .unwrap()
}

#[test]
fn barus_network_conserves_volume_and_heat_and_keeps_the_nonsymmetric_tangent() {
    let law = two_cell_barus();
    let q = [0.00001, -0.000002];
    let v = [0.0006, 0.0001];
    let mut force = [0.0; 2];
    let mut pressure = [0.0; 2];
    let report = law
        .evaluate_into(&q, &v, &mut force, &mut pressure)
        .unwrap();
    let alpha = law.pressure_viscosity_coefficient_per_pa();
    let mut outflow = [0.0; 2];
    let mut heat = 0.0;
    for edge in &law.channels {
        let h = edge.gap.reference_m - dot(&edge.gap.closure, &q);
        let g0 = edge.width_m * h.powi(3) / (12.0 * law.viscosity_pa_s().unwrap() * edge.length_m);
        let a = pressure[edge.from];
        let b = edge.to.map_or(0.0, |j| pressure[j]);
        let rate = g0 * ((-alpha * b).exp() - (-alpha * a).exp()) / alpha;
        outflow[edge.from] += rate;
        if let Some(j) = edge.to {
            outflow[j] -= rate;
        }
        heat += rate * (a - b);
    }
    for (i, cell) in law.cells.iter().enumerate() {
        assert!((outflow[i] - cell.area_m2 * dot(&cell.gap.closure, &v)).abs() < 3e-22);
        assert!((force[i] - cell.area_m2 * pressure[i]).abs() < 2e-15);
    }
    assert!((report.dissipated_power_w - heat).abs() < 3e-18);
    assert!((report.dissipated_power_w - dot(&force, &v)).abs() < 3e-18);

    // Independent centered differences vary geometry and velocity together.
    let dq = [0.000002, -0.000001];
    let dv = [0.0002, -0.0001];
    let mut tangent = [0.0; 2];
    law.tangent_into(&q, &v, &dq, &dv, &mut tangent).unwrap();
    let epsilon = 1e-5;
    let mut plus = [0.0; 2];
    let mut minus = [0.0; 2];
    let qp = [q[0] + epsilon * dq[0], q[1] + epsilon * dq[1]];
    let qm = [q[0] - epsilon * dq[0], q[1] - epsilon * dq[1]];
    let vp = [v[0] + epsilon * dv[0], v[1] + epsilon * dv[1]];
    let vm = [v[0] - epsilon * dv[0], v[1] - epsilon * dv[1]];
    law.evaluate_into(&qp, &vp, &mut plus, &mut pressure)
        .unwrap();
    law.evaluate_into(&qm, &vm, &mut minus, &mut pressure)
        .unwrap();
    for i in 0..2 {
        let finite_difference = (plus[i] - minus[i]) / (2.0 * epsilon);
        assert!((tangent[i] - finite_difference).abs() < 2e-9 * tangent[i].abs().max(0.01));
    }
    let mut first = [0.0; 2];
    let mut second = [0.0; 2];
    law.tangent_into(&q, &v, &[0.0; 2], &[1.0, 0.0], &mut first)
        .unwrap();
    law.tangent_into(&q, &v, &[0.0; 2], &[0.0, 1.0], &mut second)
        .unwrap();
    assert!((first[1] - second[0]).abs() > 0.05 * first[1].abs());
}

#[test]
fn passage_reports_close_local_volume_and_physical_work_for_heterogeneous_barus_flow() {
    let q = [0.00001, -0.000002];
    for alpha in [0.0, 2e-4] {
        let law = two_cell_barus().with_barus_viscosity(alpha).unwrap()
            .with_channel_viscosities(&[0.007, 0.02, 0.013]).unwrap();
        for sign in [-1.0, 1.0] {
            let velocity = [sign * 0.0006, sign * 0.0001];
            let mut force = [0.0; 2];
            let mut pressure = [0.0; 2];
            let mut channels = [FilmChannelReport::default(); 3];
            let report = law.evaluate_with_channels_into(
                &q, &velocity, &mut force, &mut pressure, &mut channels,
            ).unwrap();
            let mut ordinary_force = [0.0; 2];
            let mut ordinary_pressure = [0.0; 2];
            let ordinary = law.evaluate_into(
                &q, &velocity, &mut ordinary_force, &mut ordinary_pressure,
            ).unwrap();
            assert_eq!(force, ordinary_force);
            assert_eq!(pressure, ordinary_pressure);
            assert_eq!(report.dissipated_power_w, ordinary.dissipated_power_w);
            let mut outflow = [0.0; 2];
            let mut heat = 0.0;
            for (k, edge) in law.channels.iter().enumerate() {
                let h = edge.gap.reference_m - dot(&edge.gap.closure, &q);
                let g = edge.width_m * h.powi(3)
                    / (12.0 * law.channel_viscosities_pa_s()[k] * edge.length_m);
                let a = pressure[edge.from];
                let b = edge.to.map_or(0.0, |j| pressure[j]);
                // Independent integrated Poiseuille fluidity uses actual p,
                // not the pressure solver's transformed potential or factors.
                let flow = if alpha == 0.0 { g * (a-b) }
                    else { g * ((-alpha*b).exp()-(-alpha*a).exp())/alpha };
                assert!((channels[k].volume_flow_m3_s-flow).abs() < 4e-22);
                assert!((channels[k].dissipated_power_w-flow*(a-b)).abs() < 5e-18);
                assert!(channels[k].dissipated_power_w >= 0.0);
                outflow[edge.from] += channels[k].volume_flow_m3_s;
                if let Some(j) = edge.to { outflow[j] -= channels[k].volume_flow_m3_s; }
                heat += channels[k].dissipated_power_w;
            }
            for (i, cell) in law.cells.iter().enumerate() {
                assert!((outflow[i]-cell.area_m2*dot(&cell.gap.closure,&velocity)).abs() < 4e-22);
            }
            assert_eq!(heat, report.dissipated_power_w);
            assert!((heat-dot(&velocity,&force)).abs() < 5e-18);
        }
    }
}

#[test]
fn passage_reports_keep_closed_channels_zero_and_refuse_without_partial_outputs() {
    let mut law = two_cell_barus();
    law.channels[0].gap.reference_m = 0.0; // Each cell still drains to ambient.
    law.channels[0].gap.closure.fill(0.0);
    let mut force = [0.0; 2];
    let mut pressure = [0.0; 2];
    let mut channels = [FilmChannelReport::default(); 3];
    law.evaluate_with_channels_into(&[0.0;2],&[0.0001,0.0001],
        &mut force,&mut pressure,&mut channels).unwrap();
    assert_eq!(channels[0], FilmChannelReport::default());
    let sentinel = FilmChannelReport {volume_flow_m3_s:42.0,dissipated_power_w:43.0};
    for (q, velocity, count) in [
        ([0.0;2],[0.0001,0.0001],2),
        ([0.0002,0.0],[0.0;2],3),
        ([0.0;2],[1.0,0.0],3),
    ] {
        force.fill(42.0); pressure.fill(42.0); channels.fill(sentinel);
        assert!(law.evaluate_with_channels_into(&q,&velocity,
            &mut force,&mut pressure,&mut channels[..count]).is_err());
        assert_eq!(force,[42.0;2]); assert_eq!(pressure,[42.0;2]);
        assert_eq!(channels,[sentinel;3]);
    }
}

#[test]
fn zero_barus_is_bitwise_newtonian_and_pressure_barriers_refuse_atomically() {
    let original = film();
    let zero = original
        .with_barus_viscosity(2e-4)
        .unwrap()
        .with_barus_viscosity(0.0)
        .unwrap();
    let q = [0.00001, 0.0];
    let v = [0.001, 0.0];
    let mut original_f = [0.0; 2];
    let mut zero_f = [0.0; 2];
    let mut original_p = [0.0];
    let mut zero_p = [0.0];
    let a = original
        .evaluate_into(&q, &v, &mut original_f, &mut original_p)
        .unwrap();
    let b = zero
        .evaluate_into(&q, &v, &mut zero_f, &mut zero_p)
        .unwrap();
    assert_eq!(original_f, zero_f);
    assert_eq!(original_p, zero_p);
    assert_eq!(
        a.dissipated_power_w.to_bits(),
        b.dissipated_power_w.to_bits()
    );
    original
        .tangent_into(&q, &v, &[1e-6, 0.0], &[1e-4, 0.0], &mut original_f)
        .unwrap();
    zero.tangent_into(&q, &v, &[1e-6, 0.0], &[1e-4, 0.0], &mut zero_f)
        .unwrap();
    assert_eq!(original_f, zero_f);
    let law = original.with_barus_viscosity(2e-4).unwrap();
    for speed in [0.0015, 0.002] {
        // First exceeds declared pressure; second exceeds the alpha*psi=1 barrier.
        let mut force = [42.0; 2];
        let mut pressure = [42.0];
        assert!(
            law.evaluate_into(&[0.0; 2], &[speed, 0.0], &mut force, &mut pressure)
                .is_err()
        );
        assert_eq!(force, [42.0; 2]);
        assert_eq!(pressure, [42.0]);
        assert!(
            law.tangent_into(&[0.0; 2], &[speed, 0.0], &[0.0; 2], &[1.0, 0.0], &mut force)
                .is_err()
        );
        assert_eq!(force, [42.0; 2]);
    }
    for alpha in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(original.with_barus_viscosity(alpha).is_err());
    }
}

#[test]
fn barus_pressure_law_cannot_be_silently_used_by_other_film_storage_models() {
    use super::cavitation::{CavitatingFilm, CavitationLiquid, CavitationSolveConfig};
    use super::isothermal::{AmbientGas, IsothermalFilm};
    let law = film().with_barus_viscosity(2e-4).unwrap();
    let gas = AmbientGas {
        pressure_pa: 100000.0,
        temperature_k: 300.0,
        specific_gas_constant_j_kg_k: 287.0,
    };
    let liquid = CavitationLiquid {
        density_kg_m3: 1000.0,
        cavitation_pressure_pa: 95000.0,
        ambient_pressure_pa: 100000.0,
        ambient_liquid_fraction: 1.0,
    };
    let solve = CavitationSolveConfig {
        maximum_iterations: 16,
        absolute_mass_tolerance_kg: 1e-14,
        relative_mass_tolerance: 1e-10,
    };
    assert!(IsothermalFilm::at_ambient(law.clone(), &[0.0; 2], gas).is_err());
    assert!(CavitatingFilm::new(law.clone(), liquid, solve).is_err());
    let constant = law.with_barus_viscosity(0.0).unwrap();
    assert!(IsothermalFilm::at_ambient(constant.clone(), &[0.0; 2], gas).is_ok());
    assert!(CavitatingFilm::new(constant, liquid, solve).is_ok());
}

#[test]
fn g1_channel_viscosities_determine_pressure_flow_power_and_full_tangent() {
    let viscosities = [0.003, 0.017, 0.006];
    for alpha in [0.0, 2e-4] {
        let law = two_cell_barus()
            .with_barus_viscosity(alpha)
            .unwrap()
            .with_channel_viscosities(&viscosities)
            .unwrap();
        assert_eq!(law.channel_count(), 3);
        assert_eq!(law.viscosity_pa_s(), None);
        assert_eq!(law.channel_viscosities_pa_s(), viscosities);
        let q = [0.00001, -0.000002];
        let v = [0.0006, 0.0001];
        let mut force = [0.0; 2];
        let mut pressure = [0.0; 2];
        let report = law
            .evaluate_into(&q, &v, &mut force, &mut pressure)
            .unwrap();

        // Independent two-node Kirchhoff solution using actual geometry and
        // the individual supplied viscosities, not the compiled factors.
        let g: Vec<_> = law
            .channels
            .iter()
            .zip(viscosities)
            .map(|(edge, mu)| {
                let h = edge.gap.reference_m - dot(&edge.gap.closure, &q);
                edge.width_m * h.powi(3) / (12.0 * mu * edge.length_m)
            })
            .collect();
        let source = [0.0001 * v[0], 0.0002 * v[1]];
        let determinant = g[0] * g[1] + g[0] * g[2] + g[1] * g[2];
        let psi = [
            ((g[0] + g[1]) * source[0] + g[0] * source[1]) / determinant,
            (g[0] * source[0] + (g[0] + g[2]) * source[1]) / determinant,
        ];
        let exact = psi.map(|x| {
            if alpha == 0.0 {
                x
            } else {
                -(1.0 - alpha * x).ln() / alpha
            }
        });
        for i in 0..2 {
            assert!((pressure[i] - exact[i]).abs() < 2e-11);
            assert!((force[i] - law.cells[i].area_m2 * exact[i]).abs() < 2e-15);
        }
        let flux = [g[0] * (psi[0] - psi[1]), g[1] * psi[1], g[2] * psi[0]];
        assert!((flux[0] + flux[2] - source[0]).abs() < 3e-22);
        assert!((flux[1] - flux[0] - source[1]).abs() < 3e-22);
        let heat = flux[0] * (exact[0] - exact[1]) + flux[1] * exact[1] + flux[2] * exact[0];
        assert!((report.dissipated_power_w - heat).abs() < 3e-18);
        assert!((dot(&force, &v) - heat).abs() < 3e-18);

        let dq = [0.000002, -0.000001];
        let dv = [0.0002, -0.0001];
        let mut tangent = [0.0; 2];
        law.tangent_into(&q, &v, &dq, &dv, &mut tangent).unwrap();
        let epsilon = 1e-5;
        let qp: [f64; 2] = std::array::from_fn(|i| q[i] + epsilon * dq[i]);
        let qm: [f64; 2] = std::array::from_fn(|i| q[i] - epsilon * dq[i]);
        let vp: [f64; 2] = std::array::from_fn(|i| v[i] + epsilon * dv[i]);
        let vm: [f64; 2] = std::array::from_fn(|i| v[i] - epsilon * dv[i]);
        let (mut plus, mut minus) = ([0.0; 2], [0.0; 2]);
        law.evaluate_into(&qp, &vp, &mut plus, &mut pressure)
            .unwrap();
        law.evaluate_into(&qm, &vm, &mut minus, &mut pressure)
            .unwrap();
        for i in 0..2 {
            let finite_difference = (plus[i] - minus[i]) / (2.0 * epsilon);
            assert!((tangent[i] - finite_difference).abs() < 2e-9 * tangent[i].abs().max(0.01));
        }
    }
}

#[test]
fn g0_g3_uniform_channel_field_is_bitwise_original_and_invalid_updates_are_atomic() {
    let original = two_cell_barus();
    let uniform = original.with_channel_viscosities(&[0.01; 3]).unwrap();
    let restored = original
        .with_channel_viscosities(&[0.002, 0.007, 0.013])
        .unwrap()
        .with_viscosity(0.01)
        .unwrap();
    for law in [uniform, restored] {
        assert_eq!(law.viscosity_pa_s(), Some(0.01));
        assert_eq!(law.channel_viscosities_pa_s(), &[0.01; 3]);
        assert_eq!(law.factors, original.factors);
        let q = [0.00001, -0.000002];
        let v = [0.0006, 0.0001];
        let (mut a, mut b, mut pa, mut pb) = ([0.0; 2], [0.0; 2], [0.0; 2], [0.0; 2]);
        let ra = original.evaluate_into(&q, &v, &mut a, &mut pa).unwrap();
        let rb = law.evaluate_into(&q, &v, &mut b, &mut pb).unwrap();
        assert_eq!(a, b);
        assert_eq!(pa, pb);
        assert_eq!(
            ra.dissipated_power_w.to_bits(),
            rb.dissipated_power_w.to_bits()
        );
        original
            .tangent_into(&q, &v, &[1e-6, 0.0], &[1e-4, 0.0], &mut a)
            .unwrap();
        law.tangent_into(&q, &v, &[1e-6, 0.0], &[1e-4, 0.0], &mut b)
            .unwrap();
        assert_eq!(a, b);
    }
    for values in [
        vec![],
        vec![0.01; 2],
        vec![0.01; 4],
        vec![0.01, 0.0, 0.01],
        vec![0.01, -0.01, 0.01],
        vec![0.01, f64::NAN, 0.01],
        vec![0.01, f64::INFINITY, 0.01],
        vec![0.01, f64::MAX, 0.01],
    ] {
        assert!(original.with_channel_viscosities(&values).is_err());
        assert_eq!(original.channel_viscosities_pa_s(), &[0.01; 3]);
    }
}
