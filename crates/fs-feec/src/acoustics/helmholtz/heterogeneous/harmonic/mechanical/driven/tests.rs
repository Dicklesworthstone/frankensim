//! Original public FE source binding and observation integration.
use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_phs::WallPin;
use fs_rep_mesh::TetComplex;

const POINTS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const MEDIUM: Medium = Medium {
    density_kg_m3: 2.0,
    sound_speed_m_s: 3.0,
};
const FACE: [u32; 3] = [0, 1, 3];
fn uniform<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}
fn wall() -> MechanicalBoundary {
    MechanicalBoundary {
        faces: vec![FACE],
        wall: WallPin {
            surface_density: 1.0,
            stiffness_per_area: 4.0,
            resistance: 0.4,
        },
        relaxation: vec![],
    }
}
fn profile() -> WallTractionProfile {
    WallTractionProfile {
        wall: 0,
        coefficients: vec![(0, 1.0), (1, 1.0), (3, 1.0)],
    }
}

#[test]
fn original_wall_traction_drives_observed_peak_and_bandwidth_under_pressure_release() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POINTS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let wall = wall();
    let profile = profile();
    let sources = [(&profile, C64::new(0.7, -0.3))];
    let drive = MechanicalResonanceDrive {
        wall_tractions: &sources,
        ..Default::default()
    };
    // All fluid coordinates are released by two distinct nonmechanical faces;
    // the original full mechanical trace is retained on the third face.
    let release = [[0, 1, 2], [1, 2, 3]];
    let controls = ResonanceControls {
        absolute_frequency_tolerance_rad_s: 2e-8,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 180,
        measure_half_power_bandwidth: true,
    };
    let result = model
        .driven_resonance_with_coupled_mechanical_walls(
            (1.4, 2.5),
            &[],
            core::slice::from_ref(&wall),
            &release,
            drive,
            MechanicalResonanceObservation::WallVelocity {
                wall: 0,
                weights: &[1.0 / 3.0; 3],
            },
            MechanicalHarmonicLimits::default(),
            controls,
            &mut || false,
        )
        .unwrap();
    assert!((result.peak.angular_frequency_rad_s - 2.0).abs() < 3e-8);
    assert!((result.half_power.unwrap().bandwidth_rad_s - 0.4).abs() < 4e-8);
    assert!((result.observation - sources[0].1.scale(1.0 / 0.4)).abs() < 2e-7);
}

#[test]
fn coherent_pressure_and_native_volume_source_binding_are_retained() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POINTS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let foreign = uniform(&geometry);
    let source = model
        .volume_velocity_source(3, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let foreign_source = foreign
        .volume_velocity_source(3, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let volume = [(&source, C64::new(0.03, -0.02))];
    let profile = profile();
    let tractions = [(&profile, C64::new(0.7, -0.3))];
    let drive = MechanicalResonanceDrive {
        volume_sources: &volume,
        wall_tractions: &tractions,
        ..Default::default()
    };
    let points = [
        MechanicalPressureObservation {
            tetrahedron: 0,
            barycentric: [0.0, 0.0, 0.0, 1.0],
            weight: 2.0,
        },
        MechanicalPressureObservation {
            tetrahedron: 0,
            barycentric: [1.0, 0.0, 0.0, 0.0],
            weight: -3.0,
        },
    ];
    let observation = MechanicalResonanceObservation::CoherentPressure(&points);
    let omega = 1.7;
    let wall = wall();
    let result = model
        .mechanical_resonance_sample(
            omega,
            true,
            &[],
            &[wall.clone()],
            &[[0, 1, 2]],
            drive,
            observation,
            MechanicalHarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    // Independent P1 tetrahedron moments with only p3 free: volume load b3=Q/24,
    // A33=1/12-omega²/1080, W33=1/12 and uniform force integral W*1=1/6.
    let d = C64::new(4.0 - omega * omega, 0.4 * omega);
    let dd = C64::new(-2.0 * omega, 0.4);
    let a = 1.0 / 12.0 - omega * omega / 1080.0;
    let da = -omega / 540.0;
    let b = volume[0].1.scale(1.0 / 24.0);
    let denominator = d.scale(a) - C64::from_re(omega * omega / 12.0);
    let denominator_slope = d.scale(da) + dd.scale(a) - C64::from_re(omega / 6.0);
    let numerator = C64::new(0.0, omega) * b * d + tractions[0].1.scale(omega * omega / 6.0);
    let numerator_slope = C64::new(0.0, 1.0) * b * d
        + C64::new(0.0, omega) * b * dd
        + tractions[0].1.scale(omega / 3.0);
    let pressure = numerator / denominator;
    let slope = (numerator_slope - pressure * denominator_slope) / denominator;
    assert!((result.0 - pressure.scale(2.0)).abs() < 2e-10);
    assert!((result.1 - slope.scale(2.0)).abs() < 2e-10);
    let wrong = [(&foreign_source, C64::ONE)];
    let drive = MechanicalResonanceDrive {
        volume_sources: &wrong,
        ..drive
    };
    assert!(
        model
            .mechanical_resonance_sample(
                omega,
                true,
                &[],
                &[wall],
                &[[0, 1, 2]],
                drive,
                observation,
                MechanicalHarmonicLimits::default(),
                &mut || false
            )
            .is_err()
    );
}
