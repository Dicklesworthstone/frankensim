//! G0/G3: strict full-band loading, physical-reference admission, and the
//! independently sampled frequency-domain pressure/volume-flow boundary.
use super::*;
use crate::driving_point::{characteristic_line_dense, impedance_line, regional::RegionalDuct};
use fs_duct::{MouthLoad, Segment, TabulatedMouth, Termination};
use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasSpec;

const N: usize = 512;
fn fixture(rate: u32) -> (Duct, GasState, TabulatedLoad, f64) {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.005,
            length: 0.18,
        }],
    };
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.005_f64.powi(2));
    // An explicit passive resistive port. Its exact frequency independence
    // isolates table consumption from any interpolation or radiation model.
    let resistance = 0.7 * zc;
    let table = TabulatedLoad::try_new(
        vec![1.0e-4, core::f64::consts::TAU * f64::from(rate)],
        vec![C64::new(resistance, 0.0); 2],
        "resistive-terminal-control",
        0.0,
    )
    .unwrap()
    .with_mouth(TabulatedMouth {
        axial_position_m: 0.18,
        radius_m: 0.005,
        gas,
        dc_impedance: resistance,
    })
    .unwrap();
    (duct, gas, table, zc)
}

fn spectrum(mut line: DelayedFilter) -> Vec<C64> {
    let mut samples: Vec<FftC64> = (0..N)
        .map(|j| FftC64::new(line.push(if j == 0 { 1.0 } else { 0.0 }).unwrap(), 0.0))
        .collect();
    Fft::new(N).forward(&mut samples, &mut vec![FftC64::new(0.0, 0.0); N]);
    samples.iter().map(|z| C64::new(z.re, -z.im)).collect()
}

#[test]
fn full_band_characteristic_and_outlet_lines_match_frequency_domain_at_both_rates() {
    for rate in [48_000, 96_000] {
        let (duct, gas, table, zc) = fixture(rate);
        let line = characteristic_line_dense(
            &duct,
            &gas,
            MouthLoad::Tabulated(&table),
            rate,
            8192,
            zc,
            None,
            Some(N),
        )
        .unwrap();
        let actual = spectrum(line);
        let outlet = spectrum(
            impedance_line(&duct, &gas, MouthLoad::Tabulated(&table), rate, N, None).unwrap(),
        );
        for k in 1..N / 2 {
            let omega = core::f64::consts::TAU * k as f64 * f64::from(rate) / N as f64;
            let reference = fs_duct::input_impedance_tabulated(
                &duct,
                &gas,
                omega,
                LossModel::Bessel,
                &table,
                None,
            )
            .unwrap();
            let expected = fs_vfit::discretize::reflectance(reference.impedance, zc);
            assert!((actual[k] - expected).abs() < 2e-11, "rate {rate}, bin {k}");
            assert!((outlet[k] - reference.p_mouth_over_u_in).abs() < 2e-10 * zc);
        }
        let high = core::f64::consts::TAU * 0.45 * f64::from(rate);
        assert!(
            fs_duct::input_impedance(
                &duct,
                &gas,
                high,
                LossModel::Bessel,
                Termination::UnflangedOpen
            )
            .is_err(),
            "control must exceed analytic ka support"
        );
        let dc = response(
            &duct,
            &gas,
            1e-3,
            C64::new(table.mouth().unwrap().dc_impedance, 0.0),
            WallModel::Rigid,
        )
        .unwrap();
        let expected_dc = fs_vfit::discretize::reflectance(dc.impedance, zc).re;
        assert!((actual[0].re - expected_dc).abs() < 2e-12);
        assert!(
            (actual[0].re + 1.0).abs() > 0.1,
            "explicit DC must not become an ideal opening"
        );
    }
}

#[test]
fn uniform_regional_realization_is_identical_and_load_changes_are_live() {
    let (duct, gas, table, zc) = fixture(48_000);
    let regional = RegionalDuct::new(&duct, &[gas]).unwrap();
    // The scalar API takes the caller's port impedance; the regional owner
    // derives rho*c/((pi*r)*r), while the fixture divides by pi*(r*r). The two
    // areas differ by at most their own two roundings each (4u), plus one
    // division each (2u) and this check's quotient (u): 7u < 4 eps.
    // Bit identity is a claim about the realization, so share the SAME port.
    let port = regional.inlet_impedance();
    assert!((port / zc - 1.0).abs() <= 4.0 * f64::EPSILON);
    let uniform = spectrum(
        characteristic_line_dense(
            &duct,
            &gas,
            MouthLoad::Tabulated(&table),
            48_000,
            8192,
            port,
            None,
            Some(N),
        )
        .unwrap(),
    );
    let same = spectrum(
        regional
            .characteristic_line(MouthLoad::Tabulated(&table), 48_000, 8192, None, Some(N))
            .unwrap(),
    );
    assert_eq!(uniform, same);
    let changed = TabulatedLoad::try_new(
        vec![1e-4, 1e6],
        vec![C64::new(1.3 * zc, 0.0); 2],
        "different-resistive-terminal",
        0.0,
    )
    .unwrap()
    .with_mouth(TabulatedMouth {
        dc_impedance: 1.3 * zc,
        ..*table.mouth().unwrap()
    })
    .unwrap();
    let different = spectrum(
        regional
            .characteristic_line(MouthLoad::Tabulated(&changed), 48_000, 8192, None, Some(N))
            .unwrap(),
    );
    assert!((uniform[31] - different[31]).abs() > 0.05);
}

#[test]
fn missing_reference_changed_geometry_or_medium_and_incomplete_band_refuse() {
    let (duct, gas, table, zc) = fixture(48_000);
    let attempt = |table: &TabulatedLoad| {
        characteristic_line_dense(
            &duct,
            &gas,
            MouthLoad::Tabulated(table),
            48_000,
            8192,
            zc,
            None,
            Some(N),
        )
    };
    let unbound =
        TabulatedLoad::try_new(vec![1e-4, 1e6], vec![C64::new(zc, 0.0); 2], "unbound", 0.0)
            .unwrap();
    assert!(attempt(&unbound).is_err());
    for mouth in [
        TabulatedMouth {
            axial_position_m: 0.181,
            ..*table.mouth().unwrap()
        },
        TabulatedMouth {
            radius_m: 0.006,
            ..*table.mouth().unwrap()
        },
        TabulatedMouth {
            gas: GasState::try_new(&GasSpec::dry_air_ussa1976(), 303.15, 101_325.0).unwrap(),
            ..*table.mouth().unwrap()
        },
    ] {
        assert!(attempt(&table.clone().with_mouth(mouth).unwrap()).is_err());
    }
    for support in [[1000.0, 1e6], [1e-4, 1000.0]] {
        let truncated = TabulatedLoad::try_new(
            support.to_vec(),
            vec![C64::new(zc, 0.0); 2],
            "incomplete",
            0.0,
        )
        .unwrap()
        .with_mouth(*table.mouth().unwrap())
        .unwrap();
        assert!(matches!(
            attempt(&truncated),
            Err(DrivingPointError::Duct(DuctError::OutOfTable { .. }))
        ));
        assert!(
            impedance_line(
                &duct,
                &gas,
                MouthLoad::Tabulated(&truncated),
                48_000,
                N,
                None
            )
            .is_err()
        );
    }
}

#[test]
fn active_table_cannot_be_hidden_by_a_permissive_source_tolerance() {
    let (duct, gas, table, zc) = fixture(48_000);
    let active = TabulatedLoad::try_new(
        vec![1e-4, 1e6],
        vec![C64::new(-0.5 * zc, 0.0); 2],
        "deliberately-active-mutation",
        2.0,
    )
    .unwrap()
    .with_mouth(*table.mouth().unwrap())
    .unwrap();
    assert!(
        characteristic_line_dense(
            &duct,
            &gas,
            MouthLoad::Tabulated(&active),
            48_000,
            8192,
            zc,
            None,
            Some(N)
        )
        .is_err()
    );
    assert!(impedance_line(&duct, &gas, MouthLoad::Tabulated(&active), 48_000, N, None).is_err());
    // Strong bore attenuation can conceal gain at the inlet. Admission must
    // still reject the active MOUTH, independently of propagation losses.
    let long = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0001,
            length: 1.0,
        }],
    };
    let masked = active
        .with_mouth(TabulatedMouth {
            axial_position_m: 1.0,
            radius_m: 0.0001,
            ..*table.mouth().unwrap()
        })
        .unwrap();
    assert!(
        characteristic_line_dense(
            &long,
            &gas,
            MouthLoad::Tabulated(&masked),
            48_000,
            8192,
            zc,
            None,
            Some(N)
        )
        .is_err()
    );
}
