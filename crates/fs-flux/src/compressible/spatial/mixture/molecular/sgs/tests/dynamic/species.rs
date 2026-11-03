use super::*;

fn config() -> MixtureDynamicSpeciesSubgridConfig {
    MixtureDynamicSpeciesSubgridConfig {
        transport: DynamicThermalSubgridConfig {
            momentum: dynamic_config().transport.model,
            heat: DynamicScalarConfig {
                averaging: DynamicAveraging::LocalBox,
                maximum_coefficient: 100.0,
                maximum_diffusivity_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
        },
        species: DynamicComponentsConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient: 100.0,
            usage: SgsUse::Exploration,
        },
        maximum_species_diffusivity_m2_s: 10.0,
    }
}

fn species_field(formation: f64) -> MixtureEuler3d {
    let grid = CartesianGrid::new(9, 7, 5, 0.09, 0.14, 0.20).unwrap();
    field(material(formation), grid, PERIODIC, |point| {
        let [x, y, _] = std::array::from_fn(|axis| {
            std::f64::consts::TAU * point[axis] / (grid.shape[axis] as f64 * grid.spacing[axis])
        });
        let fraction = 0.5 + 0.125 * (x + y).sin();
        (
            [fraction, 1.0 - fraction, 0.0],
            [0.3 * y.sin(), 0.0, 0.0],
            1.0 + 0.1 * (x + y).cos(),
        )
    })
}

fn prepared(
    flow: &MixtureEuler3d,
    cfg: MixtureDynamicSpeciesSubgridConfig,
) -> MolecularTransport3d<ConstantTransport> {
    let source = original(flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_species_sgs(flow, &source)
            .unwrap();
    let base = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
    assert_eq!(bytes - base, 352 * flow.cell_count());
    let op =
        MolecularTransport3d::with_dynamic_species_sgs(flow, source, cfg, bytes, || false).unwrap();
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert_eq!(op.dynamic_species_subgrid_config(), Some(cfg));
    assert_eq!(op.dynamic_thermal_subgrid_config(), None);
    assert_eq!(op.dynamic_subgrid_config(), None);
    assert_eq!(op.subgrid_config(), None);
    op
}

fn fixed_schmidt(flow: &MixtureEuler3d) -> MolecularTransport3d<ConstantTransport> {
    let source = original(flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_thermal_sgs(flow, &source)
            .unwrap();
    MolecularTransport3d::with_dynamic_thermal_sgs(
        flow,
        source,
        MixtureDynamicThermalSubgridConfig {
            transport: config().transport,
            turbulent_schmidt: 0.8,
        },
        bytes,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_dynamic_species_flux_and_source_enthalpy_survive_zero_eddy_viscosity() {
    let flow = species_field(0.0);
    let mut active = prepared(&flow, config());
    let mut fixed = fixed_schmidt(&flow);
    let rate = faces(&mut active, &flow).rate;
    let base_rate = faces(&mut fixed, &flow).rate;
    assert!(rate > base_rate);
    assert!(
        active
            .sgs
            .as_ref()
            .unwrap()
            .dynamic_viscosity
            .iter()
            .all(|&mu| mu == 0.0)
    );
    let ViscosityClosure::DynamicSpecies(species) = &active.sgs.as_ref().unwrap().closure else {
        unreachable!()
    };
    let mut nonzero = 0;
    for axis in 0..3 {
        for left in 0..flow.cell_count() {
            let mut c = flow.grid.coordinates(left);
            c[axis] = (c[axis] + 1) % flow.grid.shape[axis];
            let right = flow.grid.index(c);
            let face = flow.grid.cell_face(axis, c);
            let pl = active.work.primitive[left];
            let pr = active.work.primitive[right];
            let mass_diffusion = species.face_species(left, right);
            let gradient = (flow.states()[7 * right] / pr.density
                - flow.states()[7 * left] / pl.density)
                / flow.grid.spacing[axis];
            let flux = -mass_diffusion * gradient;
            let delta: [f64; 7] = std::array::from_fn(|k| {
                active.work.flux[axis][7 * face + k] - fixed.work.flux[axis][7 * face + k]
            });
            close(delta[0], flux, 4e-14);
            close(delta[1], -flux, 4e-14);
            assert_eq!(delta[2], 0.0);
            for k in 3..6 {
                assert_eq!(delta[k], 0.0);
            }
            let at = Temperature::new(0.5 * pl.temperature + 0.5 * pr.temperature);
            let h: Vec<_> = flow
                .gas()
                .source()
                .models()
                .iter()
                .take(2)
                .map(|model| {
                    (model
                        .molar_internal_energy_change(
                            Temperature::new(flow.gas().domain().minimum_temperature),
                            at,
                        )
                        .unwrap()
                        .value()
                        + R * at.value())
                        / model.molar_mass().value()
                })
                .collect();
            close(delta[6], (h[0] - h[1]) * flux, 4e-14);
            nonzero += usize::from(flux.abs() > 1e-8);
        }
    }
    assert!(nonzero > 0);
}

#[test]
fn g3_dynamic_species_updates_conserve_mix_and_preserve_source_gauge() {
    let mut flow = species_field(0.0);
    let mut shifted = species_field(1e12);
    let mut base = species_field(0.0);
    let mut active = prepared(&flow, config());
    let mut shifted_op = prepared(&shifted, config());
    let mut fixed = fixed_schmidt(&base);
    let bytes = active.workspace_bytes() + flow.workspace_bytes();
    let controls = DissipationControl::default();
    let report = active
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    balance(&report);
    shifted_op
        .advance(&mut shifted, Time::new(0.001), controls, || false)
        .unwrap();
    fixed
        .advance(&mut base, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(flow.states(), shifted.states());
    let variance = |flow: &MixtureEuler3d| {
        flow.states()
            .chunks_exact(7)
            .map(|row| {
                let rho = sum(row[..3].iter().copied());
                rho * (row[0] / rho - 0.5).powi(2)
            })
            .sum::<f64>()
    };
    assert!(variance(&flow) < variance(&base));
    assert!(flow.states().chunks_exact(7).all(|row| row[2] == 0.0));
    assert_eq!(active.workspace_bytes() + flow.workspace_bytes(), bytes);
}

#[test]
fn g0_dynamic_species_zero_fit_clears_cache_and_preserves_molecular_bits() {
    let flow = species_field(0.0);
    let mut active = prepared(&flow, config());
    faces(&mut active, &flow);
    let uniform = field(material(0.0), flow.grid, PERIODIC, |_| {
        ([0.5, 0.5, 0.0], [0.0; 3], 1.0)
    });
    let mut base = operator(&uniform, None);
    faces(&mut active, &uniform);
    faces(&mut base, &uniform);
    let ViscosityClosure::DynamicSpecies(species) = &active.sgs.as_ref().unwrap().closure else {
        unreachable!()
    };
    assert!((0..flow.cell_count()).all(|cell| species.face_species(cell, cell) == 0.0));
    assert_eq!(active.work.flux, base.work.flux);
}

#[test]
fn g4_dynamic_species_bounds_and_late_cancellation_preserve_accepted_state() {
    let mut flow = species_field(0.0);
    let source = original(&flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_species_sgs(&flow, &source)
            .unwrap();
    assert!(
        MolecularTransport3d::with_dynamic_species_sgs(&flow, source, config(), bytes - 1, || {
            false
        },)
        .is_err()
    );
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        let mut cfg = config();
        cfg.maximum_species_diffusivity_m2_s = value;
        assert!(cfg.validate().is_err());
    }
    let mut forbidden = config();
    forbidden.species.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    let mut cfg = config();
    cfg.maximum_species_diffusivity_m2_s = 0.0;
    let mut active = prepared(&flow, cfg);
    let before = flow.states().to_vec();
    let controls = DissipationControl::default();
    assert!(
        active
            .advance(&mut flow, Time::new(0.001), controls, || false)
            .is_err()
    );
    assert_eq!(flow.states(), before);
    assert_eq!(active.last_step(), None);
    let ViscosityClosure::DynamicSpecies(species) = &mut active.sgs.as_mut().unwrap().closure
    else {
        unreachable!()
    };
    species.config = config();

    let mut reference_flow = species_field(0.0);
    let mut reference = prepared(&reference_flow, config());
    let mut polls = 0;
    let report = reference
        .advance(&mut reference_flow, Time::new(0.001), controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let loads = active.boundary_impulses.clone();
    let mut cancelled_polls = 0;
    assert_eq!(
        active.advance(&mut flow, Time::new(0.001), controls, || {
            cancelled_polls += 1;
            cancelled_polls == polls - 1
        }),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(flow.states(), before);
    assert_eq!(active.boundary_impulses, loads);
    assert_eq!(active.last_step(), None);
    let retry = active
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(retry, report);
    assert_eq!(flow.states(), reference_flow.states());
}
