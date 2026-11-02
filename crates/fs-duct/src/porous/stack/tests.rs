use super::*;
use fs_la::eigen_complex::lu_complex;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};

const AREA: f64 = 0.004;
const HZ: f64 = 1100.0;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn materials() -> [PorousFluidLaw; 3] {
    std::array::from_fn(|i| {
        let mut card = PorousCard::melamine_uf_nguyen2024();
        card.name = "synthetic heterogeneous porous stack";
        card.provenance = "numerical constitutive fixture, not measured qualification";
        card.porosity = [0.68, 0.93, 0.81][i];
        card.flow_resistivity *= [0.5, 2.0, 1.3][i];
        card.viscous_length_m *= [0.8, 1.2, 0.6][i];
        PorousFluidLaw::try_new(card, 1.81e-5).unwrap()
    })
}
fn layer(material: &PorousFluidLaw, length: f64) -> PorousLayer<'_> {
    let gas = gas();
    PorousLayer::try_new(
        material,
        &gas,
        core::f64::consts::TAU * HZ,
        length,
        AREA,
        gas.characteristic_impedance / AREA,
    )
    .unwrap()
}
fn options(n: usize) -> ComposeOptions {
    ComposeOptions {
        maximum_blocks: n,
        ..ComposeOptions::default()
    }
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-20),
        "{a:.17e} != {b:.17e}"
    );
}
fn complex_close(a: C64, b: C64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-12),
        "{a:?} != {b:?}"
    );
}
fn phase(k: C64, x: f64) -> C64 {
    C64::new((k.re * x).cos(), (k.re * x).sin()).scale((-k.im * x).exp())
}

// Independent simultaneous interface equations with all internal backward
// amplitudes anchored at each LEFT face. This intentionally uses growing
// exponentials, restricted to moderately thick reference specimens only.
struct Reference {
    outgoing: [C64; 2],
    waves: Vec<[C64; 2]>,
}
fn reference(layers: &[PorousLayer<'_>], incoming: [C64; 2]) -> Reference {
    let n = 2 + 2 * layers.len();
    let mut a = vec![C64::ZERO; n * n];
    let mut b = vec![C64::ZERO; n];
    let r = layers[0].two_port().reference_impedance;
    let impedance: Vec<_> = layers
        .iter()
        .map(|layer| C64::from_re(r * AREA) / layer.material_response().characteristic().0.conj())
        .collect();
    let travel: Vec<_> = layers
        .iter()
        .map(|layer| {
            phase(
                layer.material_response().characteristic().1.conj(),
                layer.thickness_m(),
            )
        })
        .collect();
    let one = C64::ONE;
    // Unknowns are the two exterior outgoing pressures, then F_i,B_i.
    a[0] = one;
    a[2] = -one;
    a[3] = -one;
    b[0] = -incoming[0];
    a[n] = -one;
    a[n + 2] = -impedance[0];
    a[n + 3] = impedance[0];
    b[1] = -incoming[0];
    let last = layers.len() - 1;
    let j = 2 + 2 * last;
    a[2 * n + 1] = one;
    a[2 * n + j] = -travel[last];
    a[2 * n + j + 1] = -travel[last].recip();
    b[2] = -incoming[1];
    a[3 * n + 1] = one;
    a[3 * n + j] = -impedance[last] * travel[last];
    a[3 * n + j + 1] = impedance[last] / travel[last];
    b[3] = incoming[1];
    for i in 0..last {
        let row = 4 + 2 * i;
        let col = 2 + 2 * i;
        a[row * n + col] = travel[i];
        a[row * n + col + 1] = travel[i].recip();
        a[row * n + col + 2] = -one;
        a[row * n + col + 3] = -one;
        a[(row + 1) * n + col] = impedance[i] * travel[i];
        a[(row + 1) * n + col + 1] = -impedance[i] / travel[i];
        a[(row + 1) * n + col + 2] = -impedance[i + 1];
        a[(row + 1) * n + col + 3] = impedance[i + 1];
    }
    lu_complex(&a, n).unwrap().solve(&mut b);
    Reference {
        outgoing: [b[0], b[1]],
        waves: (0..layers.len())
            .map(|i| [b[2 + 2 * i], b[3 + 2 * i]])
            .collect(),
    }
}
fn reference_field(layer: &PorousLayer<'_>, waves: [C64; 2], x: f64) -> [C64; 2] {
    let (z, k) = layer.material_response().characteristic();
    let e = phase(k.conj(), x);
    let forward = waves[0] * e;
    let backward = waves[1] / e;
    [forward + backward, (forward - backward) / z.conj()]
}
fn reference_work(layer: &PorousLayer<'_>, waves: [C64; 2]) -> [f64; 2] {
    let response = layer.material_response();
    let weights = [
        -0.5 * response.angular_frequency_rad_s() * response.dynamic_density_kg_m3().im,
        -0.5 * response.angular_frequency_rad_s() * response.dynamic_bulk_modulus_pa().recip().im,
    ];
    let panels = 1024;
    let dx = layer.thickness_m() / panels as f64;
    let mut sum = [0.0; 2];
    for i in 0..=panels {
        let [pressure, velocity] = reference_field(layer, waves, i as f64 * dx);
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum[0] += weight * weights[0] * velocity.norm_sq();
        sum[1] += weight * weights[1] * pressure.norm_sq();
    }
    sum.map(|value| value * AREA * dx / 3.0)
}

#[test]
fn g1_heterogeneous_stack_closes_interfaces_and_independent_volume_heat() {
    let materials = materials();
    let layers = [
        layer(&materials[0], 0.018),
        layer(&materials[1], 0.024),
        layer(&materials[2], 0.011),
    ];
    let stack = PorousLayerStack::try_new(&layers, options(3), || false).unwrap();
    assert_eq!(stack.layers().len(), 3);
    assert!(std::ptr::eq(stack.layers()[1].material(), &materials[1]));
    let mut powers = Vec::new();
    for incoming in [
        [C64::ONE, C64::ZERO],
        [C64::ZERO, C64::new(0.45, 0.2)],
        [C64::ONE, C64::new(0.45, 0.2)],
    ] {
        let excitation = stack.prepare_excitation(incoming, || false).unwrap();
        let independent = reference(&layers, incoming);
        for (actual, expected) in excitation
            .boundary_power()
            .outgoing_pressure_pa
            .into_iter()
            .zip(independent.outgoing)
        {
            complex_close(actual, expected, 2e-12);
        }
        for i in 0..3 {
            let local = &excitation.layer_excitations()[i];
            for depth in [0.0, 0.37 * layers[i].thickness_m(), layers[i].thickness_m()] {
                let actual = local.at_depth(depth).unwrap();
                let expected = reference_field(&layers[i], independent.waves[i], depth);
                complex_close(actual.pressure_pa, expected[0], 3e-12);
                complex_close(actual.superficial_velocity_m_s, expected[1], 3e-12);
            }
            let work = local
                .dissipation_between(0.0, layers[i].thickness_m())
                .unwrap();
            let expected = reference_work(&layers[i], independent.waves[i]);
            close(work.viscous_power_w, expected[0], 3e-11);
            close(work.thermal_power_w, expected[1], 3e-11);
            if i > 0 {
                let left = excitation.layer_excitations()[i - 1]
                    .at_depth(layers[i - 1].thickness_m())
                    .unwrap();
                let right = local.at_depth(0.0).unwrap();
                complex_close(left.pressure_pa, right.pressure_pa, 3e-12);
                complex_close(left.volume_flow_m3_s, right.volume_flow_m3_s, 3e-12);
            }
        }
        close(
            excitation.dissipation().absorbed_power_w,
            excitation.boundary_power().absorbed_power_w,
            2e-12,
        );
        let whole = excitation
            .dissipation_between(0.0, stack.thickness_m())
            .unwrap();
        assert_eq!(whole, excitation.dissipation());
        let selected = excitation.at_depth(stack.layer_depths_m()[1]).unwrap();
        assert_eq!(selected.layer_index, 1);
        assert_eq!(selected.local_depth_m, 0.0);
        powers.push(whole.absorbed_power_w);
    }
    // Both coherent sources are combined before squaring, so the shared drive
    // is not a sum of independent one-sided heat fractions.
    assert!((powers[2] - powers[0] - powers[1]).abs() > 1e-3 * powers[2]);
}

#[test]
fn g3_one_layer_and_homogeneous_splits_preserve_fields_and_cell_integrals() {
    let materials = materials();
    let original = layer(&materials[0], 0.08);
    let one = [original.clone()];
    let split = [
        layer(&materials[0], 0.013),
        layer(&materials[0], 0.017),
        layer(&materials[0], 0.05),
    ];
    for layers in [&one[..], &split[..]] {
        let stack = PorousLayerStack::try_new(layers, options(3), || false).unwrap();
        for incoming in [
            [C64::ONE, C64::ZERO],
            [C64::ZERO, C64::ONE],
            [C64::new(0.3, -0.2), C64::new(-0.4, 0.7)],
        ] {
            let native = original.prepare_excitation(incoming).unwrap();
            let excitation = stack.prepare_excitation(incoming, || false).unwrap();
            for depth in [0.0, 0.005, 0.013, 0.02, 0.03, 0.045, 0.08] {
                let actual = excitation.at_depth(depth).unwrap().field;
                let expected = native.at_depth(depth).unwrap();
                complex_close(actual.pressure_pa, expected.pressure_pa, 3e-12);
                complex_close(
                    actual.superficial_velocity_m_s,
                    expected.superficial_velocity_m_s,
                    3e-12,
                );
            }
            for (start, end) in [(0.0, 0.021), (0.02, 0.055), (0.0, 0.08), (0.03, 0.03)] {
                let actual = excitation.dissipation_between(start, end).unwrap();
                let expected = native.dissipation_between(start, end).unwrap();
                close(actual.viscous_power_w, expected.viscous_power_w, 4e-12);
                close(actual.thermal_power_w, expected.thermal_power_w, 4e-12);
                close(actual.absorbed_power_w, expected.absorbed_power_w, 4e-12);
            }
            for (actual, expected) in excitation
                .boundary_power()
                .outgoing_pressure_pa
                .into_iter()
                .zip(native.boundary_power().outgoing_pressure_pa)
            {
                complex_close(actual, expected, 3e-12);
            }
        }
    }
}

#[test]
fn g1_opaque_two_sided_stack_and_physical_reversal_keep_face_heat() {
    let materials = materials();
    let layers = [
        layer(&materials[0], 1000.0),
        layer(&materials[1], 0.03),
        layer(&materials[2], 1000.0),
    ];
    let reversed: Vec<_> = layers.iter().rev().cloned().collect();
    let stack = PorousLayerStack::try_new(&layers, options(3), || false).unwrap();
    let reverse = PorousLayerStack::try_new(&reversed, options(3), || false).unwrap();
    let incoming = [C64::ONE, C64::new(0.4, -0.7)];
    let excitation = stack.prepare_excitation(incoming, || false).unwrap();
    let backward = reverse
        .prepare_excitation([incoming[1], incoming[0]], || false)
        .unwrap();
    assert_eq!(stack.two_port().scattering[1], C64::ZERO);
    assert_eq!(stack.two_port().scattering[2], C64::ZERO);
    assert_eq!(
        excitation.layer_excitations()[1]
            .dissipation_between(0.0, 0.03)
            .unwrap()
            .absorbed_power_w,
        0.0
    );
    assert_eq!(
        excitation.at_depth(500.0).unwrap().field.pressure_pa,
        C64::ZERO
    );
    assert_eq!(
        excitation
            .at_depth(1500.0)
            .unwrap()
            .field
            .superficial_velocity_m_s,
        C64::ZERO
    );
    for i in 0..3 {
        let actual = excitation.layer_excitations()[i]
            .dissipation_between(0.0, layers[i].thickness_m())
            .unwrap();
        let mirrored = backward.layer_excitations()[2 - i]
            .dissipation_between(0.0, layers[i].thickness_m())
            .unwrap();
        close(actual.absorbed_power_w, mirrored.absorbed_power_w, 3e-12);
        if i != 1 {
            let local_drive = if i == 0 {
                [incoming[0], C64::ZERO]
            } else {
                [C64::ZERO, incoming[1]]
            };
            let local = layers[i].prepare_excitation(local_drive).unwrap();
            close(
                actual.absorbed_power_w,
                local
                    .dissipation_between(0.0, layers[i].thickness_m())
                    .unwrap()
                    .absorbed_power_w,
                3e-12,
            );
        }
    }
    close(
        excitation.dissipation().absorbed_power_w,
        excitation.boundary_power().absorbed_power_w,
        3e-12,
    );
}

#[test]
fn g0_stack_common_domains_work_bounds_cancellation_and_depth_refusals() {
    let materials = materials();
    let gas = gas();
    let original = layer(&materials[0], 0.02);
    assert!(matches!(
        PorousLayerStack::try_new(&[], options(2), || false),
        Err(DuctError::EmptyDuct)
    ));
    for (frequency, area, reference) in [
        (HZ + 1.0, AREA, gas.characteristic_impedance / AREA),
        (HZ, 2.0 * AREA, gas.characteristic_impedance / AREA),
        (HZ, AREA, 1.1 * gas.characteristic_impedance / AREA),
    ] {
        let changed = PorousLayer::try_new(
            &materials[1],
            &gas,
            core::f64::consts::TAU * frequency,
            0.03,
            area,
            reference,
        )
        .unwrap();
        let layers = [original.clone(), changed];
        assert!(PorousLayerStack::try_new(&layers, options(2), || false).is_err());
    }
    let layers = [original.clone(), original];
    assert!(PorousLayerStack::try_new(&layers, options(1), || false).is_err());
    assert!(
        PorousLayerStack::try_new(
            &layers,
            ComposeOptions {
                maximum_condition: f64::NAN,
                ..options(2)
            },
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        PorousLayerStack::try_new(&layers, options(2), || {
            polls += 1;
            polls > 3
        }),
        Err(DuctError::Cancelled { .. })
    ));
    let stack = PorousLayerStack::try_new(&layers, options(2), || false).unwrap();
    assert!(matches!(
        stack.prepare_excitation([C64::ONE; 2], || true),
        Err(DuctError::Cancelled { .. })
    ));
    assert!(
        stack
            .prepare_excitation([C64::new(f64::NAN, 0.0), C64::ZERO], || false)
            .is_err()
    );
    assert!(
        stack
            .prepare_excitation([C64::from_re(f64::MAX); 2], || false)
            .is_err()
    );
    let zero = stack.prepare_excitation([C64::ZERO; 2], || false).unwrap();
    assert_eq!(zero.dissipation().absorbed_power_w, 0.0);
    assert_eq!(zero.power_residual_w(), 0.0);
    for depth in [-1.0, f64::NAN, f64::INFINITY, stack.thickness_m() + 1.0] {
        assert!(zero.at_depth(depth).is_err());
        assert!(zero.dissipation_between(0.0, depth).is_err());
    }
    assert!(zero.dissipation_between(0.02, 0.01).is_err());
    assert_eq!(
        zero.dissipation_between(0.02, 0.02)
            .unwrap()
            .absorbed_power_w,
        0.0
    );
}
