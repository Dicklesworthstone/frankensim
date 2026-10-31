use super::*;
use crate::volumes::transport::{CellPalette, CellPaletteBudget};
use std::sync::Arc;

fn prepared(cx: &Cx<'_>, family: Family, temperature: f64, host: f64) -> MaterialParticleSpectrum {
    binding(&card(family), temperature, host)
        .visible_spectrum(
            cx,
            &[380.0, 550.0, 780.0],
            MiePopulationBudget::default(),
            0.003,
            2048,
        )
        .unwrap()
}

fn palette(
    cx: &Cx<'_>,
    materials: Vec<MaterialParticleSpectrum>,
    indices: &[u32],
) -> CellPalette<MaterialParticleSpectrum> {
    let max_models = materials.len();
    CellPalette::try_new(
        cx,
        materials,
        Arc::from(indices),
        CellPaletteBudget {
            max_models,
            max_cells: indices.len(),
        },
    )
    .unwrap()
}

#[test]
fn g1_two_material_cells_use_native_extinction_phase_and_one_shared_map() {
    with_cx(|cx, _| {
        let materials = palette(
            cx,
            vec![
                prepared(cx, Family::Cauchy, 300.0, 1.0),
                prepared(cx, Family::Conductor, 400.0, 1.0),
            ],
            &[0, 1],
        );
        let original = materials.clone();
        let density = [2e-5, 7e-5];
        let field =
            MaterialParticleCellGrid::try_new(cx, grid(&density), materials, None, 1e-4, budget())
                .unwrap();
        assert_eq!(field.materials().values(), original.values());
        let absorption = field
            .medium()
            .absorption_grid()
            .absorption_palette()
            .unwrap();
        let scattering = field.medium().scattering_palette().unwrap();
        assert!(Arc::ptr_eq(
            original.cell_indices(),
            absorption.cell_indices()
        ));
        assert!(Arc::ptr_eq(
            original.cell_indices(),
            scattering.cell_indices()
        ));
        assert!(Arc::ptr_eq(
            original.cell_indices(),
            field.materials().cell_indices()
        ));
        assert!(!field.medium().has_thermal_emission());
        let wavelengths = [400.0, 550.0, 760.0];
        let actual = field
            .medium()
            .transmittance(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                22.0,
                &wavelengths,
            )
            .unwrap();
        for (got, wavelength) in actual.into_iter().zip(wavelengths) {
            let mut tau = 0.0;
            for (cell, &density) in density.iter().enumerate() {
                let spectrum = &original.at_cell(cell).unwrap().spectrum;
                tau += 10.0
                    * density
                    * (spectrum.absorption.at(wavelength).unwrap()
                        + spectrum.scattering.at(wavelength).unwrap());
                for cosine in [-1.0, 0.0, 0.7] {
                    close(
                        field
                            .medium()
                            .phase_at_cell(cell)
                            .unwrap()
                            .eval_at(wavelength, cosine)
                            .unwrap(),
                        spectrum.phase.eval_at(wavelength, cosine).unwrap(),
                        2e-14,
                    );
                }
            }
            close(got, (-tau).exp(), 3e-14);
        }
        close(
            field.maximum_volume_fraction(),
            density[1] / original.values()[1].density_kg_m3,
            1e-15,
        );
    });
}

#[test]
fn g1_indexed_material_temperature_and_bulk_density_control_lte_and_dilution() {
    with_cx(|cx, _| {
        let materials = palette(
            cx,
            vec![
                prepared(cx, Family::Cauchy, 300.0, 1.0),
                prepared(cx, Family::Cauchy, 400.0, 1.0),
            ],
            &[0, 1],
        );
        let density = [2e-5, 7e-5];
        let temperatures = [300.0, 400.0];
        let field = MaterialParticleCellGrid::try_new(
            cx,
            grid(&density),
            materials.clone(),
            Some(&temperatures),
            1e-4,
            budget(),
        )
        .unwrap();
        assert!(field.medium().has_thermal_emission());
        let wavelength = 760.0;
        let transport = field
            .medium()
            .absorption_grid()
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                22.0,
                &[wavelength],
            )
            .unwrap();
        let tau = [0, 1].map(|cell| {
            10.0 * density[cell]
                * materials
                    .at_cell(cell)
                    .unwrap()
                    .spectrum
                    .absorption
                    .at(wavelength)
                    .unwrap()
        });
        let planck = temperatures.map(|temperature| {
            crate::thermal_surface::planck_radiance_w_m2_sr_nm(temperature, wavelength).unwrap()
        });
        let expected =
            planck[0] * (-(-tau[0]).exp_m1()) + (-tau[0]).exp() * planck[1] * (-(-tau[1]).exp_m1());
        close(transport.emission_w_m2_sr_nm[0], expected, 5e-13);
        assert!(matches!(
            MaterialParticleCellGrid::try_new(
                cx,
                grid(&density),
                materials.clone(),
                Some(&[300.0, 300.0]),
                1e-4,
                budget()
            ),
            Err(ParticleScatteringError::TemperatureMismatch { index: 1 })
        ));
        // 1000 kg/m³ in the first cell permits this density; the second
        // accepted material is 900 kg/m³ and must reject the same dispersed mass.
        assert!(matches!(
            MaterialParticleCellGrid::try_new(
                cx,
                grid(&[0.0095, 0.0095]),
                materials,
                None,
                1e-5,
                budget()
            ),
            Err(ParticleScatteringError::VolumeFraction { index: 1 })
        ));
    });
}

#[test]
fn g0_indexed_particle_grid_refuses_map_host_budget_and_cancelled_publication() {
    with_cx(|cx, gate| {
        let spectrum = prepared(cx, Family::Cauchy, 300.0, 1.0);
        let density = [2e-5, 7e-5];
        let short = palette(cx, vec![spectrum.clone()], &[0]);
        assert!(matches!(
            MaterialParticleCellGrid::try_new(cx, grid(&density), short, None, 1e-4, budget()),
            Err(ParticleScatteringError::Transport(
                VolumeTransportError::Palette
            ))
        ));
        let foreign_host = palette(cx, vec![prepared(cx, Family::Cauchy, 300.0, 1.33)], &[0, 0]);
        assert!(matches!(
            MaterialParticleCellGrid::try_new(
                cx,
                grid(&density),
                foreign_host,
                None,
                1e-4,
                budget()
            ),
            Err(ParticleScatteringError::HostIndex)
        ));
        let materials = palette(cx, vec![spectrum], &[0, 0]);
        assert!(matches!(
            MaterialParticleCellGrid::try_new(
                cx,
                grid(&density),
                materials.clone(),
                None,
                1e-4,
                AbsorbingGridBudget {
                    max_grid_cells: 1,
                    ..budget()
                }
            ),
            Err(ParticleScatteringError::Transport(
                VolumeTransportError::Budget
            ))
        ));
        gate.request();
        assert!(matches!(
            MaterialParticleCellGrid::try_new(cx, grid(&density), materials, None, 1e-4, budget()),
            Err(ParticleScatteringError::Transport(
                VolumeTransportError::Cancelled
            ))
        ));
    });
}
