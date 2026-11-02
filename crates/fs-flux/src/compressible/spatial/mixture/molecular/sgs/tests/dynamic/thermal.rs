use super::*;
use fs_lbm::sgs3::dynamic::scalar::PeriodicDynamicScalar;

fn config() -> MixtureDynamicThermalSubgridConfig {
    MixtureDynamicThermalSubgridConfig {
        transport: DynamicThermalSubgridConfig {
            momentum: dynamic_config().transport.model,
            heat: DynamicScalarConfig {
                averaging: DynamicAveraging::LocalBox,
                maximum_coefficient: 100.0,
                maximum_diffusivity_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
        },
        turbulent_schmidt: 0.8,
    }
}

fn heat_field(formation: f64) -> MixtureEuler3d {
    heat_field_on_grid(
        formation,
        CartesianGrid::new(9, 9, 9, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0).unwrap(),
    )
}

fn heat_field_on_grid(formation: f64, grid: CartesianGrid) -> MixtureEuler3d {
    field(material(formation), grid, PERIODIC, |point| {
        let [x, y, _] = std::array::from_fn(|axis| {
            std::f64::consts::TAU * point[axis] / (grid.shape[axis] as f64 * grid.spacing[axis])
        });
        // Exact-zero momentum fit: its Lxx is orthogonal to Mxy. The
        // oblique temperature field has a nonzero scalar Germano fit.
        (
            [0.6, 0.4, 0.0],
            [0.3 * y.sin(), 0.0, 0.0],
            1.0 + 0.1 * (x + y).sin(),
        )
    })
}

fn prepared(
    flow: &MixtureEuler3d,
    cfg: MixtureDynamicThermalSubgridConfig,
) -> MolecularTransport3d<ConstantTransport> {
    let source = original(flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_thermal_sgs(flow, &source)
            .unwrap();
    let base = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
    assert_eq!(bytes - base, 264 * flow.cell_count());
    let op =
        MolecularTransport3d::with_dynamic_thermal_sgs(flow, source, cfg, bytes, || false).unwrap();
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert_eq!(op.dynamic_thermal_subgrid_config(), Some(cfg));
    assert_eq!(op.dynamic_subgrid_config(), None);
    assert_eq!(op.subgrid_config(), None);
    op
}

#[test]
fn g1_dynamic_mixture_heat_survives_zero_eddy_viscosity_and_uses_original_cp() {
    let flow = heat_field(0.0);
    check_heat_faces(&flow);
}

fn check_heat_faces(flow: &MixtureEuler3d) -> (MolecularTransport3d<ConstantTransport>, f64) {
    let cfg = config();
    let mut op = prepared(flow, cfg);
    let mut molecular = operator(flow, None);
    let dynamic_rate = faces(&mut op, flow).rate;
    let molecular_rate = faces(&mut molecular, flow).rate;
    assert!(dynamic_rate > molecular_rate);
    assert!(
        op.sgs
            .as_ref()
            .unwrap()
            .dynamic_viscosity
            .iter()
            .all(|&mu| mu == 0.0)
    );
    let rho: Vec<_> = op.work.primitive.iter().map(|p| p.density).collect();
    let velocity: Vec<_> = op.work.primitive.iter().map(|p| p.velocity).collect();
    let temperature: Vec<_> = op.work.primitive.iter().map(|p| p.temperature).collect();
    let mut core =
        PeriodicDynamicSgs::new_cartesian(flow.grid.shape, flow.grid.spacing, usize::MAX).unwrap();
    let mut heat = PeriodicDynamicScalar::new(flow.grid.shape, usize::MAX).unwrap();
    let (momentum, scalar) = core
        .prepare_with_scalar(
            cfg.transport.momentum,
            cfg.transport.heat,
            &rho,
            &velocity,
            &temperature,
            &mut heat,
            || false,
        )
        .unwrap();
    assert_eq!(momentum.maximum_eddy_viscosity, 0.0);
    assert!(scalar.maximum_diffusivity > 0.0);
    let mut active = 0;
    for axis in 0..3 {
        for left in 0..flow.cell_count() {
            let mut coordinates = flow.grid.coordinates(left);
            coordinates[axis] = (coordinates[axis] + 1) % flow.grid.shape[axis];
            let right = flow.grid.index(coordinates);
            let face = flow.grid.cell_face(axis, coordinates);
            let pl = op.work.primitive[left];
            let pr = op.work.primitive[right];
            let at = Temperature::new(0.5 * pl.temperature + 0.5 * pr.temperature);
            let cp: f64 = flow
                .gas()
                .source()
                .models()
                .iter()
                .zip([0.6, 0.4, 0.0])
                .map(|(model, y)| {
                    y / model.molar_mass().value()
                        * model.evaluate_properties(at).unwrap().cp().value()
                })
                .sum();
            let weight = 0.5 * rho[left] * heat.diffusivities()[left]
                + 0.5 * rho[right] * heat.diffusivities()[right];
            let expected =
                -weight * cp * (pr.temperature - pl.temperature) / flow.grid.spacing[axis];
            for component in 0..6 {
                assert_eq!(
                    op.work.flux[axis][7 * face + component].to_bits(),
                    molecular.work.flux[axis][7 * face + component].to_bits()
                );
            }
            let actual = op.work.flux[axis][7 * face + 6] - molecular.work.flux[axis][7 * face + 6];
            close(actual, expected, 4e-14);
            active += usize::from(expected.abs() > 1e-8);
        }
    }
    assert!(active > 0);
    (op, dynamic_rate)
}

#[test]
fn g3_dynamic_mixture_rectangular_heat_uses_directional_stability_and_conserves() {
    let grid = CartesianGrid::new(9, 7, 5, 0.09, 0.14, 0.20).unwrap();
    let thin_grid = CartesianGrid::new(9, 7, 5, 0.09, 0.14, 0.05).unwrap();
    let mut flow = heat_field_on_grid(0.0, grid);
    let thin = heat_field_on_grid(0.0, thin_grid);
    let (mut op, rate) = check_heat_faces(&flow);
    let (_, thin_rate) = check_heat_faces(&thin);

    // The resolved state is independent of z. Changing only hz leaves the
    // dimensional fitted heat coefficient and molecular diffusivity unchanged;
    // the owner must still include the narrower axis in its diffusion bound.
    let inverse_square_sum =
        |grid: CartesianGrid| grid.spacing.iter().map(|h| h.recip().powi(2)).sum::<f64>();
    assert!(thin_rate > rate);
    close(
        thin_rate / rate,
        inverse_square_sum(thin_grid) / inverse_square_sum(grid),
        2e-12,
    );

    let mut molecular_flow = heat_field_on_grid(0.0, grid);
    let mut molecular = operator(&molecular_flow, None);
    let bytes = op.workspace_bytes() + flow.workspace_bytes();
    let controls = DissipationControl::default();
    let report = op
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    balance(&report);
    molecular
        .advance(&mut molecular_flow, Time::new(0.001), controls, || false)
        .unwrap();
    assert_ne!(flow.states(), molecular_flow.states());
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert!(flow.states().chunks_exact(7).all(|row| row[2] == 0.0));
}

#[test]
fn g3_dynamic_mixture_thermal_stage_conserves_and_preserves_caloric_reference() {
    let mut flow = dynamic_field(0.0);
    let mut shifted = dynamic_field(1e12);
    let mut base = dynamic_field(0.0);
    let mut op = prepared(&flow, config());
    let mut shifted_op = prepared(&shifted, config());
    let mut fixed_prandtl = dynamic_operator(&base, dynamic_config());
    let controls = DissipationControl::default();
    let report = op
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    balance(&report);
    shifted_op
        .advance(&mut shifted, Time::new(0.001), controls, || false)
        .unwrap();
    fixed_prandtl
        .advance(&mut base, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(flow.states(), shifted.states());
    assert_ne!(flow.states(), base.states());
    for row in flow.states().chunks_exact(7) {
        assert_eq!(row[2], 0.0);
    }
}

#[test]
fn g0_dynamic_mixture_thermal_zero_fit_clears_heat_cache_and_retains_molecular_bits() {
    let active = heat_field(0.0);
    let mut op = prepared(&active, config());
    faces(&mut op, &active);
    let uniform = field(material(0.0), active.grid, PERIODIC, |_| {
        ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
    });
    let mut base = operator(&uniform, None);
    faces(&mut op, &uniform);
    faces(&mut base, &uniform);
    let ViscosityClosure::DynamicThermal(thermal) = &op.sgs.as_ref().unwrap().closure else {
        unreachable!()
    };
    assert!((0..uniform.cell_count()).all(|i| thermal.face_heat(i, i).0 == 0.0));
    assert_eq!(op.work.flux, base.work.flux);
    let mut actual = uniform;
    let mut expected = field(material(0.0), actual.grid, PERIODIC, |_| {
        ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
    });
    let controls = DissipationControl::default();
    let a = op
        .advance(&mut actual, Time::new(0.001), controls, || false)
        .unwrap();
    let b = base
        .advance(&mut expected, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(actual.states(), expected.states());
}

#[test]
fn g4_dynamic_mixture_heat_budget_ratio_and_cancellation_preserve_accepted_state() {
    let mut flow = heat_field(0.0);
    let source = original(&flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_thermal_sgs(&flow, &source)
            .unwrap();
    assert!(
        MolecularTransport3d::with_dynamic_thermal_sgs(&flow, source, config(), bytes - 1, || {
            false
        },)
        .is_err()
    );
    let mut forbidden = config();
    forbidden.transport.heat.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    let mut op = prepared(&flow, config());
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
                polls == 1000
            }
        ),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(flow.states(), state);
    assert_eq!(op.boundary_impulses, loads);
    assert_eq!(op.last_step(), Some(&report));
    let ViscosityClosure::DynamicThermal(thermal) = &mut op.sgs.as_mut().unwrap().closure else {
        unreachable!()
    };
    thermal
        .config
        .transport
        .heat
        .maximum_diffusivity_to_molecular = 0.0;
    assert_eq!(
        op.advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false
        ),
        Err(DissipationError::Invalid(
            "mixture subgrid heat/molecular conductivity ratio exceeded"
        ))
    );
    assert_eq!(flow.states(), state);
    assert_eq!(op.last_step(), Some(&report));
}
