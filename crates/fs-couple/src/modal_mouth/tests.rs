use super::*;
use fs_bem::radiation_bake::lathe_profile;
use fs_bem::reciprocal_radiation;
use fs_material::gas::GasSpec;
use fs_math::c64::C64;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn budget() -> ModalMouthBemBudget {
    ModalMouthBemBudget {
        maximum_panels: 10_000,
        maximum_trace_entries: 100_000,
        maximum_matrix_entries: 1000,
        maximum_projection_terms: 1_000_000,
        maximum_basis_gram_defect: 0.2,
    }
}
fn fixture(around: usize, radial: usize) -> (SpherePanels, CircularMouthPatch) {
    let radius = 0.02;
    let profile: Vec<_> = (0..=2 * radial)
        .map(|i| (i as f64 * radius / radial as f64, radius))
        .collect();
    let bell = lathe_profile(&profile, around).unwrap();
    let patch = CircularMouthPatch {
        centre_m: [2.0 * radius, 0.0, 0.0],
        outward_normal: [1.0, 0.0, 0.0],
        radius_m: radius,
        panels: (bell.mouth_start..bell.triangles.len()).collect(),
    };
    (SpherePanels::from_triangles(bell.triangles).unwrap(), patch)
}
fn relative(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-25),
        "{a:?} != {b:?}"
    );
}

#[test]
fn signed_higher_mode_traces_retain_zero_flow_limit_and_exact_discrete_work() {
    let gate = CancelGate::new();
    let mut defects = Vec::new();
    let mut net_higher = Vec::new();
    for (around, radial) in [(12, 3), (24, 6), (48, 12)] {
        let (surface, patch) = fixture(around, radial);
        let binding =
            BemModalMouth::new(&surface, patch, air(), ModeCount(2), budget(), &gate).unwrap();
        defects.push(binding.basis_gram_defect());
        let area = core::f64::consts::PI * binding.patch().radius_m.powi(2);
        for &panel in &binding.patch().panels {
            assert!(
                (binding.velocity_trace(0).unwrap()[panel] * area - 1.0).abs() < 4.0 * f64::EPSILON
            );
        }
        let higher = binding.velocity_trace(1).unwrap();
        assert!(higher.iter().any(|v| *v < 0.0));
        assert!(higher.iter().any(|v| *v > 0.0));
        net_higher.push(
            higher
                .iter()
                .zip(surface.areas())
                .map(|(v, a)| v * a)
                .sum::<f64>()
                .abs(),
        );
        // An independently authored complex pressure and modal velocity bank.
        // The same spatial trace must be transposed for work, even when net
        // higher-mode volume flow tends to zero.
        let flow = [C64::new(0.31, -0.08), C64::new(-0.73, 0.19)];
        let mut pressure = [C64::ZERO; 2];
        let mut surface_power = 0.0;
        for i in 0..surface.areas().len() {
            let p = C64::new(
                1.1 + 19.0 * surface.centroids()[i][1],
                0.7 - 3.0 * surface.centroids()[i][2],
            );
            let mut v = C64::ZERO;
            for m in 0..2 {
                let trace = binding.velocity_trace(m).unwrap()[i];
                v = v + flow[m].scale(trace);
                pressure[m] = pressure[m] + p.scale(surface.areas()[i] * trace);
            }
            surface_power += 0.5 * surface.areas()[i] * (p * v.conj()).re;
        }
        let port_power = 0.5 * (pressure[0] * flow[0].conj() + pressure[1] * flow[1].conj()).re;
        assert!((surface_power - port_power).abs() < 1e-12);
    }
    assert!(
        defects[1] < defects[0] * 0.4 && defects[2] < defects[1] * 0.4,
        "{defects:?}"
    );
    assert!(
        net_higher[1] < net_higher[0] * 0.4 && net_higher[2] < net_higher[1] * 0.4,
        "{net_higher:?}"
    );
}

#[test]
fn complete_native_galerkin_matrix_enters_duct_without_net_flow_rescaling() {
    let (surface, patch) = fixture(12, 3);
    let gate = CancelGate::new();
    let binding =
        BemModalMouth::new(&surface, patch, air(), ModeCount(2), budget(), &gate).unwrap();
    let route = ApertureRoute::Galerkin {
        budget: reciprocal_radiation::Budget {
            maximum_unknowns: 1000,
            maximum_cache_bytes: 64 << 20,
            maximum_work: 4_000_000_000,
            angular_refinement: 0,
        },
    };
    let omegas = [1100.0, 2300.0];
    let baked = binding.bake(&omegas, &route, &gate).unwrap();
    assert_eq!(baked.loads().len(), 2);
    assert_eq!(baked.basis_gram_defect(), binding.basis_gram_defect());
    assert_eq!(baked.cap_area_ratio(), binding.cap_area_ratio());
    let gain = &baked.source().volume_gains;
    assert!((gain[0] - binding.cap_area_ratio()).abs() < 1e-14);
    assert!(gain[1].abs() < 0.08);
    assert!(binding.velocity_trace(1).unwrap().iter().any(|v| *v < 0.0));
    let duct = fs_duct::Duct {
        segments: vec![fs_duct::Segment::Cylinder {
            radius: binding.patch().radius_m,
            length: 0.04,
        }],
    };
    for (load, row) in baked.loads().iter().zip(&baked.source().rows) {
        assert_eq!(row.route, "galerkin");
        assert_eq!(load.omega_rad_s(), row.omega);
        assert_eq!(load.impedance(), row.impedance.as_slice());
        assert!(row.impedance[1].abs() > 1e-8 * row.impedance[0].abs());
        relative(row.impedance[1], row.impedance[2], 1e-10);
        for u in [
            [C64::ONE, C64::ZERO],
            [C64::ZERO, C64::ONE],
            [C64::new(0.2, -0.6), C64::new(0.9, 0.3)],
        ] {
            let p = [
                row.impedance[0] * u[0] + row.impedance[1] * u[1],
                row.impedance[2] * u[0] + row.impedance[3] * u[1],
            ];
            let power = 0.5 * (u[0].conj() * p[0] + u[1].conj() * p[1]).re;
            assert!(power >= -1e-12 * row.impedance[0].abs());
        }
        // Falsifier: converting this signed basis to independent plane-flow
        // coordinates changes the full matrix and is not the product path.
        let wrong = baked
            .source()
            .acoustic(
                baked
                    .loads()
                    .iter()
                    .position(|v| v.omega_rad_s() == row.omega)
                    .unwrap(),
            )
            .unwrap();
        assert!((wrong[3] - row.impedance[3]).abs() > row.impedance[3].abs());
        // Complete physical path: unchanged exterior -> loaded bore -> the
        // actual modal mouth flow -> independent exterior boundary work.
        let response = fs_duct::modal::outlet::mouth_response_matrix(
            &duct,
            &air(),
            row.omega,
            fs_duct::LossModel::Lossless,
            load,
            2,
            1,
        )
        .unwrap();
        let u = [response.flow_map[0], response.flow_map[2]]; // U_in=[1,0]
        let p = [
            row.impedance[0] * u[0] + row.impedance[1] * u[1],
            row.impedance[2] * u[0] + row.impedance[3] * u[1],
        ];
        let outgoing = 0.5 * (u[0].conj() * p[0] + u[1].conj() * p[1]).re;
        let supplied = 0.5 * response.input.plane_impedance.re;
        assert!((outgoing - supplied).abs() < 1e-8 * outgoing.abs().max(supplied.abs()));
        let diagonal = ModalMouthImpedance::new(
            row.omega,
            binding.patch().radius_m,
            air(),
            ModeCount(2),
            vec![row.impedance[0], C64::ZERO, C64::ZERO, row.impedance[3]],
        )
        .unwrap();
        let stripped = fs_duct::modal::mm_input_impedance_matrix(
            &duct,
            &air(),
            row.omega,
            fs_duct::LossModel::Lossless,
            &diagonal,
            2,
            1,
        )
        .unwrap();
        assert!(
            (stripped.plane_impedance - response.input.plane_impedance).abs()
                > 1e-5 * response.input.plane_impedance.abs()
        );
    }
}

#[test]
fn invalid_geometry_unresolved_basis_budget_and_cancellation_refuse() {
    let (surface, patch) = fixture(16, 3);
    let gate = CancelGate::new();
    for malformed in 0..5 {
        let mut bad_patch = patch.clone();
        match malformed {
            0 => bad_patch.panels.push(bad_patch.panels[0]),
            1 => bad_patch.outward_normal[0] = -1.0,
            2 => bad_patch.centre_m[0] += 1e-5,
            3 => bad_patch.radius_m *= 0.9,
            _ => {
                bad_patch.panels.pop();
            }
        }
        let mut options = budget();
        // Missing cap coverage is a physical basis discrepancy, not an index
        // error. A tight declared Gram allowance must expose it.
        if malformed == 4 {
            options.maximum_basis_gram_defect = 1e-6;
        }
        assert!(
            BemModalMouth::new(&surface, bad_patch, air(), ModeCount(2), options, &gate).is_err()
        );
    }
    for kind in 0..4 {
        let mut options = budget();
        match kind {
            0 => options.maximum_panels = 1,
            1 => options.maximum_trace_entries = 1,
            2 => options.maximum_projection_terms = 1,
            _ => options.maximum_basis_gram_defect = 1e-6,
        }
        assert!(
            BemModalMouth::new(&surface, patch.clone(), air(), ModeCount(2), options, &gate)
                .is_err()
        );
    }
    let mut options = budget();
    options.maximum_matrix_entries = 1;
    let bound =
        BemModalMouth::new(&surface, patch.clone(), air(), ModeCount(2), options, &gate).unwrap();
    let route = ApertureRoute::Collocation { iterative: None };
    for frequencies in [&[][..], &[100.0, 50.0][..], &[f64::NAN][..], &[100.0][..]] {
        assert!(matches!(
            bound.bake(frequencies, &route, &gate),
            Err(ModalMouthBemError::Invalid(_))
        ));
    }
    gate.request();
    assert!(matches!(
        bound.bake(&[100.0], &route, &gate),
        Err(ModalMouthBemError::Cancelled)
    ));
    assert!(matches!(
        BemModalMouth::new(&surface, patch, air(), ModeCount(2), budget(), &gate),
        Err(ModalMouthBemError::Cancelled)
    ));
}
