//! G0/G3: unsupported transverse motion cannot disappear during composition.
use super::*;
use fs_scenario::{BendingRelaxationBranch, Listener, PronyBending};

fn assembly(prony: bool) -> AcousticAssembly {
    AcousticAssembly {
        ambient: AmbientGas::sea_level(),
        string: Some(PrestressedString {
            length_m: 0.5,
            tension_n: 64.0,
            lin_density_kg_m: 1.0,
            axial_stiffness_n: 0.0,
            width_m: 0.001,
            n_modes: 1,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 1e-3,
            kelvin_voigt_bending: None,
            relaxation_bending: prony.then_some(PronyBending {
                branches: vec![BendingRelaxationBranch {
                    relaxing_stiffness_n_m2: 1e-3,
                    relaxation_time_s: 0.1,
                }],
                omega_band_rad_s: (0.0, 1e5),
                material_state_identity: None,
                source_properties: None,
            }),
            polarization_detune: 0.003,
            moving_end: true,
            tension_loss_factor: None,
            bending_loss_factor: None,
        }),
        duct: None,
        pluck: Some(Pluck {
            station_frac: 0.5,
            height_m: 1e-5,
        }),
        bow: None,
        blow: None,
        reed: None,
        soundboard: Some(RadiatingPlate {
            area_m2: 0.01,
            mass_kg: 0.2,
            frequency_hz: 10.0,
            damping_ratio: 0.02,
        }),
        body_modes: vec![],
        plate: None,
        cavity: None,
        obstacles: vec![],
        contact_texture: None,
        listener: Listener { distance_m: 1.0 },
        sample_rate_hz: 48_000,
        duration_s: 32.0 / 48_000.0,
    }
}

#[test]
fn moving_end_cannot_silently_omit_a_second_polarization() {
    let expected = AcousticRealizeError::InvalidDescription {
        what: "moving-end second polarization requires explicit transverse attachment geometry and is not supported by this assembly",
    };
    for prony in [false, true] {
        for body in [false, true] {
            let mut a = assembly(prony);
            if !body {
                a.soundboard = None;
            }
            assert_eq!(realize_assembly(&a).unwrap_err(), expected);
            // Refusal belongs to the assembly topology, before material/gas
            // realization can mask it or allocate a partially realized body.
            a.ambient.pressure_pa = -1.0;
            assert_eq!(realize_assembly(&a).unwrap_err(), expected);
        }
    }
}

#[test]
fn supported_single_plane_and_fixed_end_polarizations_still_render() {
    for prony in [false, true] {
        for body in [false, true] {
            for (moving_end, detune) in [(true, 0.0), (false, 0.003)] {
                let mut a = assembly(prony);
                if !body {
                    a.soundboard = None;
                }
                let string = a.string.as_mut().unwrap();
                string.moving_end = moving_end;
                string.polarization_detune = detune;
                let pressure = realize_assembly(&a).unwrap().pressure_pa;
                assert_eq!(pressure.len(), 32);
                assert!(pressure.iter().all(|p| p.is_finite()));
                assert!(pressure.iter().any(|p| *p != 0.0));
            }
        }
    }
}
