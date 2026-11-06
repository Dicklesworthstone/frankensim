use super::*;
use fs_lbm::sgs3::dynamic::PeriodicDynamicSgs;

mod channel;
mod species;
mod thermal;

fn dynamic_config() -> MixtureDynamicSubgridConfig {
    MixtureDynamicSubgridConfig {
        transport: DynamicSubgridConfig {
            model: DynamicSgsConfig {
                averaging: DynamicAveraging::LocalBox,
                maximum_coefficient_squared: 100.0,
                maximum_eddy_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
            turbulent_prandtl: 0.7,
        },
        turbulent_schmidt: 0.8,
    }
}

fn dynamic_field(formation: f64) -> MixtureEuler3d {
    dynamic_field_on_grid(
        formation,
        CartesianGrid::new(7, 7, 7, 1.0 / 7.0, 1.0 / 7.0, 1.0 / 7.0).unwrap(),
    )
}

fn dynamic_field_on_grid(formation: f64, grid: CartesianGrid) -> MixtureEuler3d {
    field(material(formation), grid, PERIODIC, |point| {
        let [x, y, z] = std::array::from_fn(|axis| {
            std::f64::consts::TAU * point[axis] / (grid.shape[axis] as f64 * grid.spacing[axis])
        });
        let rho = 1.2 + 0.2 * x.sin() + 0.1 * y.cos();
        let fraction = 0.5 + 0.15 * (x + y).sin();
        (
            [rho * fraction, rho * (1.0 - fraction), 0.0],
            [
                0.3 * x.sin() + 0.2 * (y + z).sin() + 0.1 * (2.0 * x - y).cos(),
                0.2 * y.cos() + 0.3 * (z + x).sin() + 0.1 * (x - 2.0 * z).sin(),
                0.2 * z.sin() + 0.1 * (x + y).cos() + 0.3 * (2.0 * y - z).sin(),
            ],
            1.0 + 0.05 * z.cos(),
        )
    })
}

fn dynamic_operator(
    flow: &MixtureEuler3d,
    config: MixtureDynamicSubgridConfig,
) -> MolecularTransport3d<ConstantTransport> {
    let source = original(flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_sgs(flow, &source).unwrap();
    let op = MolecularTransport3d::with_dynamic_sgs(flow, source, config, bytes, || false).unwrap();
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert_eq!(op.dynamic_subgrid_config(), Some(config));
    assert_eq!(op.subgrid_config(), None);
    op
}

#[test]
fn g0_dynamic_mixture_uses_actual_favre_density_and_conservative_species_faces() {
    let flow = dynamic_field(0.0);
    check_dynamic_species_faces(&flow);
}

fn check_dynamic_species_faces(flow: &MixtureEuler3d) -> MolecularTransport3d<ConstantTransport> {
    let cfg = dynamic_config();
    let mut op = dynamic_operator(flow, cfg);
    let mut molecular = operator(flow, None);
    let modeled_rate = faces(&mut op, flow).rate;
    let molecular_rate = faces(&mut molecular, flow).rate;
    assert!(modeled_rate > molecular_rate);
    let rho: Vec<_> = op.work.primitive.iter().map(|p| p.density).collect();
    let velocity: Vec<_> = op.work.primitive.iter().map(|p| p.velocity).collect();
    let mut reference =
        PeriodicDynamicSgs::new_cartesian(flow.grid.shape, flow.grid.spacing, usize::MAX).unwrap();
    let fit = reference
        .prepare(cfg.transport.model, &rho, &velocity, || false)
        .unwrap();
    assert!(fit.maximum_eddy_viscosity > 1e-5);
    let modeled = &op.sgs.as_ref().unwrap().dynamic_coefficient_squared;
    for (&coefficient, fitted) in modeled.iter().zip(reference.coefficients_squared()) {
        assert_eq!(coefficient, fitted);
    }
    reference
        .prepare(
            cfg.transport.model,
            &vec![1.0; rho.len()],
            &velocity,
            || false,
        )
        .unwrap();
    assert!(
        modeled
            .iter()
            .zip(reference.coefficients_squared())
            .any(|(&coefficient, reynolds)| (coefficient - reynolds).abs() > 1e-7),
        "the variable-density fixture must expose a Reynolds/Favre substitution"
    );

    let mut active_faces = 0;
    for axis in 0..3 {
        for left in 0..flow.cell_count() {
            let mut c = flow.grid.coordinates(left);
            c[axis] = (c[axis] + 1) % flow.grid.shape[axis];
            let right = flow.grid.index(c);
            let face = flow.grid.cell_face(axis, c);
            let mu = dynamic_face_viscosity_reference(
                &op,
                [modeled[left], modeled[right]],
                reference.filter_width(),
                left,
                right,
                axis,
            );
            let pl = op.work.primitive[left];
            let pr = op.work.primitive[right];
            let yl = flow.states()[7 * left] / pl.density;
            let yr = flow.states()[7 * right] / pr.density;
            let j = -(mu / cfg.turbulent_schmidt) * ((yr - yl) / flow.grid.spacing[axis]);
            let delta =
                |k| op.work.flux[axis][7 * face + k] - molecular.work.flux[axis][7 * face + k];
            close(delta(0), j, 2e-14);
            close(delta(1), -j, 2e-14);
            assert_eq!(delta(2), 0.0);
            active_faces += usize::from(j.abs() > 1e-8);
        }
    }
    assert!(active_faces > 0);
    op
}

#[test]
fn g3_dynamic_mixture_rectangular_cells_conserve_species_energy_and_source_gauge() {
    let grid = CartesianGrid::new(5, 7, 9, 0.16, 0.09, 0.12).unwrap();
    let mut flow = dynamic_field_on_grid(0.0, grid);
    let mut shifted = dynamic_field_on_grid(1e12, grid);
    let mut laminar = dynamic_field_on_grid(0.0, grid);
    let mut op = check_dynamic_species_faces(&flow);
    let mut shifted_op = dynamic_operator(&shifted, dynamic_config());
    let mut molecular = operator(&laminar, None);
    let bytes = op.workspace_bytes() + flow.workspace_bytes();
    let controls = DissipationControl::default();
    let report = op
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    balance(&report);
    shifted_op
        .advance(&mut shifted, Time::new(0.001), controls, || false)
        .unwrap();
    molecular
        .advance(&mut laminar, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(flow.states(), shifted.states());
    assert_ne!(flow.states(), laminar.states());
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert!(flow.states().chunks_exact(7).all(|row| row[2] == 0.0));
}

#[test]
fn g0_dynamic_mixture_zero_fit_clears_previous_stress_and_keeps_molecular_bits() {
    let flow = dynamic_field(0.0);
    let mut op = dynamic_operator(&flow, dynamic_config());
    faces(&mut op, &flow);
    assert!(
        op.sgs
            .as_ref()
            .unwrap()
            .dynamic_coefficient_squared
            .iter()
            .any(|&mu| mu > 0.0)
    );
    let uniform = field(material(0.0), flow.grid, PERIODIC, |_| {
        ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
    });
    let mut molecular = operator(&uniform, None);
    faces(&mut op, &uniform);
    faces(&mut molecular, &uniform);
    assert!(
        op.sgs
            .as_ref()
            .unwrap()
            .dynamic_coefficient_squared
            .iter()
            .all(|&mu| mu == 0.0)
    );
    assert_eq!(op.work.flux, molecular.work.flux);

    let mut actual = uniform;
    let mut expected = field(material(0.0), actual.grid, PERIODIC, |_| {
        ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
    });
    let controls = DissipationControl::default();
    let report = op
        .advance(&mut actual, Time::new(0.001), controls, || false)
        .unwrap();
    let reference = molecular
        .advance(&mut expected, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(report, reference);
    assert_eq!(actual.states(), expected.states());
}

#[test]
fn g3_dynamic_mixture_stage_updates_conserve_species_energy_and_source_gauge() {
    let cfg = dynamic_config();
    let mut flow = dynamic_field(0.0);
    let mut shifted = dynamic_field(1e12);
    let mut laminar = dynamic_field(0.0);
    let mut op = dynamic_operator(&flow, cfg);
    let mut shifted_op = dynamic_operator(&shifted, cfg);
    let mut laminar_op = operator(&laminar, None);
    let bytes = op.workspace_bytes() + flow.workspace_bytes();
    let controls = DissipationControl::default();
    let report = op
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    balance(&report);
    shifted_op
        .advance(&mut shifted, Time::new(0.001), controls, || false)
        .unwrap();
    laminar_op
        .advance(&mut laminar, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(flow.states(), shifted.states());
    assert_ne!(flow.states(), laminar.states());
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
    for row in flow.states().chunks_exact(7) {
        assert_eq!(row[2], 0.0);
    }

    let mut small_sc = cfg;
    small_sc.turbulent_schmidt = 0.0002;
    let mut ordinary = dynamic_operator(&flow, cfg);
    let mut mixing = dynamic_operator(&flow, small_sc);
    assert!(faces(&mut mixing, &flow).rate > 5.0 * faces(&mut ordinary, &flow).rate);
}

#[test]
fn g4_dynamic_mixture_budgets_topology_and_cancel_refuse_without_publication() {
    let mut flow = dynamic_field(0.0);
    let cfg = dynamic_config();
    let source = original(&flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_sgs(&flow, &source).unwrap();
    assert!(
        MolecularTransport3d::with_dynamic_sgs(&flow, source, cfg, bytes - 1, || false).is_err()
    );
    for (grid, boundary) in [
        (CartesianGrid::new(7, 7, 7, 0.1, 0.1, 0.1).unwrap(), WALLS),
        (
            CartesianGrid::new(7, 2, 7, 0.1, 0.1, 0.1).unwrap(),
            PERIODIC,
        ),
    ] {
        let bad = field(material(0.0), grid, boundary, |_| {
            ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
        });
        assert!(
            MolecularTransport3d::with_dynamic_sgs(&bad, original(&bad), cfg, usize::MAX, || false)
                .is_err()
        );
    }
    let mut forbidden = cfg;
    forbidden.transport.model.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut config = cfg;
        config.turbulent_schmidt = bad;
        assert!(config.validate().is_err());
    }

    let mut op = dynamic_operator(&flow, cfg);
    let report = op
        .advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let state = flow.states().to_vec();
    let loads = op.boundary_impulses.clone();
    let mut polls = 0;
    assert_eq!(
        op.advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || {
                polls += 1;
                // Source admission consumes several polls per cell; this reaches
                // the actual interval rather than cancelling only at constructor.
                polls == 1000
            }
        ),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(flow.states(), state);
    assert_eq!(op.boundary_impulses, loads);
    assert_eq!(op.last_step(), Some(&report));
    let ViscosityClosure::Dynamic(dynamic) = &mut op.sgs.as_mut().unwrap().closure else {
        unreachable!()
    };
    dynamic.config.transport.model.maximum_eddy_to_molecular = 0.0;
    assert!(
        op.advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false
        )
        .is_err()
    );
    assert_eq!(flow.states(), state);
    assert_eq!(op.last_step(), Some(&report));
}

#[test]
fn dynamic_face_strain_drives_mixture_stress_and_species_at_centered_nulls() {
    let n = 6;
    let at = |x: usize, y: usize, z: usize| (z * n + y) * n + x;
    let mut random = 1u32;
    let mut velocity = vec![[0.0; 3]; n * n * n];
    for u in &mut velocity {
        for component in u {
            random = random.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            *component = 0.03 * (2.0 * (random >> 8) as f64 / 16_777_216.0 - 1.0);
        }
    }
    for x in 0..4 {
        velocity[at(x, 1, 1)] = [0.0, if x % 2 == 0 { -0.01 } else { 0.01 }, 0.0];
    }
    for x in [1, 2] {
        for [y, z] in [[0, 1], [2, 1], [1, 0], [1, 2]] {
            velocity[at(x, y, z)] = [0.0; 3];
        }
    }
    let mut flow = field(
        material(0.0),
        CartesianGrid::new(n, n, n, 1.0, 1.0, 1.0).unwrap(),
        PERIODIC,
        |p| {
            let cell = at(p[0] as usize, p[1] as usize, p[2] as usize);
            let rho = 1.0 + 0.1 * (cell % 7) as f64 / 6.0;
            let fraction = 0.4 + 0.05 * (cell % 5) as f64;
            ([rho * fraction, rho * (1.0 - fraction), 0.0], velocity[cell], 1.5)
        },
    );
    let mut cfg = dynamic_config();
    cfg.transport.model.averaging = DynamicAveraging::Domain;
    let mut active = dynamic_operator(&flow, cfg);
    let mut base = operator(&flow, None);
    let base_rate = faces(&mut base, &flow).rate;
    let active_rate = faces(&mut active, &flow).rate;
    assert!(active_rate > base_rate);
    let left = at(1, 1, 1);
    let right = at(2, 1, 1);
    for cell in [left, right] {
        assert!(active.work.gradient[cell].iter().flatten().all(|g| g.abs() < 1e-15));
    }
    let coefficients = &active.sgs.as_ref().unwrap().dynamic_coefficient_squared;
    // Independent physical-unit Germano reference for the fixed LCG field.
    let expected_coefficient = 0.079_783_248_957_537_02;
    close(coefficients[left], expected_coefficient, 2e-13);
    close(coefficients[right], expected_coefficient, 2e-13);
    let mu = 1.025 * expected_coefficient * 0.02;
    let face = 7 * flow.grid.cell_face(0, [2, 1, 1]);
    let delta = |k| active.work.flux[0][face + k] - base.work.flux[0][face + k];
    close(delta(4), 0.02 * mu, 2e-13);
    let species = -(mu / cfg.turbulent_schmidt) * 0.05;
    close(delta(0), species, 2e-13);
    close(delta(1), -species, 2e-13);
    assert_eq!(delta(2), 0.0);
    let report = active
        .advance(&mut flow, Time::new(0.001), DissipationControl::default(), || false)
        .unwrap();
    balance(&report);
    assert!(flow.states().chunks_exact(7).all(|row| row[2] == 0.0));
}
