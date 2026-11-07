//! Exact conservative remap of one uniform-density planar material column.
//!
//! Original material coordinate x maps to s*x, s=m0/m1. The pressure stage's
//! affine velocity u(x)=(s-1)*x/dt is carried with that SAME material. This leaf
//! neither solves pressure nor replaces native collision/enthalpy advancement.

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct SlabCell {
    pub index: usize,
    pub fill_fraction: f64,
    pub mass_kg: f64,
    pub momentum_kg_m_s: f64,
    /// Exact integral of the carried affine velocity's kinetic storage.
    /// The native cell-mean kinetic energy can be lower; that difference is
    /// a remap defect, never deposited as material heat.
    pub carried_kinetic_energy_j: f64,
}

#[derive(Debug, PartialEq)]
pub(super) struct SlabRemap {
    pub cells: Vec<SlabCell>,
    pub height_m: f64,
    pub source_mass_kg: f64,
    pub source_momentum_kg_m_s: f64,
    pub source_kinetic_energy_j: f64,
}

/// Integrate the inverse affine overlap of each occupied destination cell.
///
/// Depth is measured outward from the anchored wall. The caller maps this
/// coordinate/index and momentum sign into its chosen Cartesian orientation.
/// Original depth includes the actual terminal-cell fraction and may be less
/// than one cell. Original and endpoint masses are kg per FULL Cartesian cell.
/// `capacity_cells` is the caller's admitted allocation/work ceiling.
///
/// Every destination receives only original material; contraction admits no
/// gas-side donor. A partial final cell retains its actual overlap fraction.
/// Uniform specific enthalpy therefore travels with the same returned masses.
/// Source totals are analytic original-material integrals, independently of
/// destination summation. The caller audits native initialization roundoff,
/// pressure/remap work and projection defects against explicit physical limits.
///
/// Returns no partial result on invalid/unrepresentable input, overflow of the
/// available depth, allocation failure, or requested cancellation. Cancellation
/// is polled at entry, every 256 destination cells, and before publication.
pub(super) fn remap(
    original_depth_cells: f64,
    capacity_cells: usize,
    dx_m: f64,
    dt_s: f64,
    original_cell_mass_kg: f64,
    endpoint_cell_mass_kg: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Option<SlabRemap> {
    if cancelled() {
        return None;
    }
    if capacity_cells == 0
        || !original_depth_cells.is_finite()
        || original_depth_cells <= 0.0
        || original_depth_cells > capacity_cells as f64
        || [dx_m, dt_s, original_cell_mass_kg, endpoint_cell_mass_kg]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return None;
    }
    // A partial original column can have a representable endpoint even when
    // the density ratio alone overflows. Combine the original depth before
    // restoring the binary exponent; do not form that ratio independently.
    let depth = product_ratio(
        [original_depth_cells, original_cell_mass_kg, 1.0, 1.0],
        [endpoint_cell_mass_kg, 1.0],
    );
    if !depth.is_finite() || depth <= 0.0 || depth > capacity_cells as f64 {
        return None;
    }
    let count = depth.ceil() as usize;
    if count == 0 || count > capacity_cells {
        return None;
    }
    let initial_height = product_ratio([original_depth_cells, dx_m, 1.0, 1.0], [1.0, 1.0]);
    let height_m = product_ratio([depth, dx_m, 1.0, 1.0], [1.0, 1.0]);
    // Preserve changes that disappear in (rounded stretch - 1), and avoid an
    // intermediate inverse tiny dt or overflowing density/time product.
    let edge_velocity = product_ratio(
        [
            original_cell_mass_kg - endpoint_cell_mass_kg,
            original_depth_cells,
            dx_m,
            1.0,
        ],
        [endpoint_cell_mass_kg, dt_s],
    );
    let source_mass_kg = product_ratio(
        [original_depth_cells, original_cell_mass_kg, 1.0, 1.0],
        [1.0, 1.0],
    );
    let source_momentum_kg_m_s =
        product_ratio([source_mass_kg, edge_velocity, 0.5, 1.0], [1.0, 1.0]);
    let source_kinetic_energy_j = product_ratio(
        [source_mass_kg, edge_velocity, edge_velocity, 1.0 / 6.0],
        [1.0, 1.0],
    );
    if [initial_height, height_m, source_mass_kg]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || [
            edge_velocity,
            source_momentum_kg_m_s,
            source_kinetic_energy_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
    {
        return None;
    }
    let mut cells = Vec::new();
    cells.try_reserve_exact(count).ok()?;
    for index in 0..count {
        if index % 256 == 0 && cancelled() {
            return None;
        }
        let fill_fraction = (depth - index as f64).min(1.0);
        if !fill_fraction.is_finite() || fill_fraction <= 0.0 {
            return None;
        }
        let mass_kg = product_ratio([endpoint_cell_mass_kg, fill_fraction, 1.0, 1.0], [1.0, 1.0]);
        // Coordinates normalized by the COMPLETE target material depth are
        // also the original material-coordinate fractions under x -> s*x.
        let width = fill_fraction / depth;
        let midpoint = index as f64 / depth + 0.5 * width;
        let momentum_kg_m_s = product_ratio([mass_kg, edge_velocity, midpoint, 1.0], [1.0, 1.0]);
        let mean_square = midpoint * midpoint + width * width / 12.0;
        let carried_kinetic_energy_j = product_ratio(
            [mass_kg, edge_velocity, edge_velocity, 0.5 * mean_square],
            [1.0, 1.0],
        );
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !momentum_kg_m_s.is_finite()
            || !carried_kinetic_energy_j.is_finite()
        {
            return None;
        }
        cells.push(SlabCell {
            index,
            fill_fraction,
            mass_kg,
            momentum_kg_m_s,
            carried_kinetic_energy_j,
        });
    }
    if cancelled() {
        return None;
    }
    Some(SlabRemap {
        cells,
        height_m,
        source_mass_kg,
        source_momentum_kg_m_s,
        source_kinetic_energy_j,
    })
}

// Finite factors and positive finite denominators. Binary scaling prevents
// density*volume, squared speed and inverse tiny dt from overflowing before a
// representable final mass/momentum/energy is formed.
pub(super) fn product_ratio(factors: [f64; 4], denominators: [f64; 2]) -> f64 {
    fn parts(mut value: f64) -> (f64, i32) {
        let mut shift = 0;
        if value < f64::MIN_POSITIVE {
            value *= 18_014_398_509_481_984.0;
            shift = -54;
        }
        let bits = value.to_bits();
        (
            f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52)),
            ((bits >> 52) & 0x7ff) as i32 - 1023 + shift,
        )
    }
    let mut sign = 1.0;
    let mut fraction = 1.0;
    let mut exponent = 0_i32;
    for value in factors {
        if value == 0.0 {
            return 0.0;
        }
        if value < 0.0 {
            sign = -sign;
        }
        let (part, power) = parts(value.abs());
        fraction *= part;
        exponent += power;
        if fraction >= 2.0 {
            fraction *= 0.5;
            exponent += 1;
        }
    }
    for denominator in denominators {
        let (part, power) = parts(denominator);
        fraction /= part;
        exponent -= power;
        if fraction < 1.0 {
            fraction *= 2.0;
            exponent -= 1;
        }
    }
    let value = if exponent > 1023 {
        f64::INFINITY
    } else if exponent >= -1022 {
        fraction * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent < -1075 {
        0.0
    } else {
        let scale = f64::from_bits(((exponent + 2045) as u64) << 52);
        (fraction * scale) * f64::MIN_POSITIVE
    };
    sign * value
}

#[cfg(test)]
mod tests {
    use super::*;

    fn remap(
        original_cells: usize,
        capacity_cells: usize,
        dx_m: f64,
        dt_s: f64,
        original_cell_mass_kg: f64,
        endpoint_cell_mass_kg: f64,
    ) -> Option<SlabRemap> {
        super::remap(
            original_cells as f64,
            capacity_cells,
            dx_m,
            dt_s,
            original_cell_mass_kg,
            endpoint_cell_mass_kg,
            || false,
        )
    }

    #[test]
    fn cancellation_never_publishes_a_partial_column() {
        assert!(super::remap(1.0, 2, 1.0, 1.0, 1.0, 1.0, || true).is_none());
        let mut checks = 0;
        assert!(
            super::remap(1.0, 512, 1.0, 1.0, 512.0, 1.0, || {
                checks += 1;
                checks >= 3
            })
            .is_none()
        );
        assert_eq!(checks, 3);
    }

    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 4e-14 * a.abs().max(b.abs()).max(1e-300),
            "{a:e} != {b:e}"
        );
    }
    fn totals(r: &SlabRemap) {
        close(r.cells.iter().map(|c| c.mass_kg).sum(), r.source_mass_kg);
        close(
            r.cells.iter().map(|c| c.momentum_kg_m_s).sum(),
            r.source_momentum_kg_m_s,
        );
        close(
            r.cells.iter().map(|c| c.carried_kinetic_energy_j).sum(),
            r.source_kinetic_energy_j,
        );
    }

    #[test]
    fn expansion_preserves_material_momentum_and_exact_kinetic_integrals() {
        let r = remap(2, 8, 0.5, 0.25, 2.0, 1.0).unwrap();
        assert_eq!(r.cells.len(), 4);
        close(r.height_m, 2.0);
        close(r.source_mass_kg, 4.0);
        close(r.source_momentum_kg_m_s, 8.0);
        close(r.source_kinetic_energy_j, 32.0 / 3.0);
        for (i, c) in r.cells.iter().enumerate() {
            let u = i as f64 + 0.5;
            close(c.fill_fraction, 1.0);
            close(c.mass_kg, 1.0);
            close(c.momentum_kg_m_s, u);
            close(c.carried_kinetic_energy_j, 0.5 * (u * u + 1.0 / 12.0));
        }
        totals(&r);
        let native: f64 = r
            .cells
            .iter()
            .map(|c| 0.5 * c.momentum_kg_m_s.powi(2) / c.mass_kg)
            .sum();
        close(r.source_kinetic_energy_j - native, 1.0 / 6.0);
    }

    #[test]
    fn contraction_transports_original_liquid_without_a_gas_mass_donor() {
        let r = remap(4, 8, 0.5, 0.25, 1.0, 2.0).unwrap();
        assert_eq!(r.cells.len(), 2);
        close(r.height_m, 1.0);
        close(r.source_mass_kg, 4.0);
        close(r.source_momentum_kg_m_s, -8.0);
        close(r.source_kinetic_energy_j, 32.0 / 3.0);
        for (i, c) in r.cells.iter().enumerate() {
            let u = -(2.0 * i as f64 + 1.0);
            close(c.mass_kg, 2.0);
            close(c.momentum_kg_m_s, 2.0 * u);
            close(c.carried_kinetic_energy_j, u * u + 1.0 / 3.0);
        }
        totals(&r);
    }

    #[test]
    fn partial_endpoint_and_arbitrary_enthalpy_datum_retain_inventory() {
        let r = remap(2, 5, 1.0, 1.0, 3.0, 4.0).unwrap();
        assert_eq!(r.cells.len(), 2);
        close(r.cells[0].mass_kg, 4.0);
        close(r.cells[1].fill_fraction, 0.5);
        close(r.cells[1].mass_kg, 2.0);
        close(r.cells[0].momentum_kg_m_s, -2.0 / 3.0);
        close(r.cells[1].momentum_kg_m_s, -5.0 / 6.0);
        close(r.source_kinetic_energy_j, 0.25);
        for h in [0.0, -10_000.0, 1e12] {
            close(
                r.cells.iter().map(|c| c.mass_kg * h).sum(),
                r.source_mass_kg * h,
            );
        }
        totals(&r);
    }

    #[test]
    fn equal_density_and_small_resolved_expansion_have_no_hidden_floor() {
        let same = remap(3, 3, 0.1, 0.01, 0.5, 0.5).unwrap();
        assert_eq!(same.cells.len(), 3);
        for c in same.cells {
            assert_eq!(c.momentum_kg_m_s, 0.0);
            assert_eq!(c.carried_kinetic_energy_j, 0.0);
        }
        let delta = 2.0_f64.powi(-40);
        let tiny = remap(1, 2, 1.0, 1.0, 1.0 + delta, 1.0).unwrap();
        assert_eq!(tiny.cells.len(), 2);
        close(tiny.cells[1].fill_fraction, delta);
        close(tiny.source_momentum_kg_m_s, (1.0 + delta) * delta / 2.0);
        assert!(tiny.source_kinetic_energy_j > 0.0);
        totals(&tiny);
    }

    #[test]
    fn density_length_and_time_rescaling_follow_physical_dimensions() {
        let a = remap(3, 9, 0.2, 0.1, 0.7, 0.4).unwrap();
        let b = remap(3, 9, 0.6, 0.2, 7.0, 4.0).unwrap();
        assert_eq!(a.cells.len(), b.cells.len());
        close(b.height_m, 3.0 * a.height_m);
        for (a, b) in a.cells.iter().zip(&b.cells) {
            close(b.fill_fraction, a.fill_fraction);
            close(b.mass_kg, 10.0 * a.mass_kg);
            close(b.momentum_kg_m_s, 15.0 * a.momentum_kg_m_s);
            close(
                b.carried_kinetic_energy_j,
                22.5 * a.carried_kinetic_energy_j,
            );
        }
        totals(&a);
        totals(&b);
    }

    #[test]
    fn extreme_representable_mass_speed_products_and_tiny_time_remain_finite() {
        for (m, dx) in [(1e-200, 1e100), (1e200, 1e-100)] {
            let r = remap(1, 2, dx, 1.0, 2.0 * m, m).unwrap();
            close(r.source_kinetic_energy_j, 1.0 / 3.0);
            totals(&r);
        }
        let tiny = remap(1, 2, 1e-310, 1e-320, 2.0, 1.0).unwrap();
        close(tiny.source_momentum_kg_m_s, 1e-310 / 1e-320);
        totals(&tiny);
    }

    #[test]
    fn invalid_inputs_capacity_and_unrepresentable_totals_refuse() {
        assert!(remap(0, 1, 1.0, 1.0, 1.0, 1.0).is_none());
        assert!(remap(2, 1, 1.0, 1.0, 1.0, 1.0).is_none());
        assert!(remap(2, 2, 1.0, 1.0, 2.0, 1.0).is_none());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(remap(1, 2, bad, 1.0, 1.0, 1.0).is_none());
            assert!(remap(1, 2, 1.0, bad, 1.0, 1.0).is_none());
            assert!(remap(1, 2, 1.0, 1.0, bad, 1.0).is_none());
            assert!(remap(1, 2, 1.0, 1.0, 1.0, bad).is_none());
        }
        assert!(remap(2, 2, 1.0, 1.0, f64::MAX, f64::MAX).is_none());
    }
    #[test]
    fn partial_source_preserves_exact_affine_mass_momentum_and_kinetic_storage() {
        // H0=1.25 m, M=7.5 kg, s=3/4, surface velocity=-1.25 m/s.
        let r = super::remap(2.5, 8, 0.5, 0.25, 3.0, 4.0, || false).unwrap();
        close(r.height_m, 0.9375);
        close(r.source_mass_kg, 7.5);
        close(r.source_momentum_kg_m_s, -4.6875);
        close(r.source_kinetic_energy_j, 1.953125);
        assert_eq!(r.cells.len(), 2);
        close(r.cells[0].fill_fraction, 1.0);
        close(r.cells[1].fill_fraction, 0.875);
        close(r.cells[0].momentum_kg_m_s, -4.0 / 3.0);
        close(r.cells[1].momentum_kg_m_s, -3.3541666666666665);
        totals(&r);
    }

    #[test]
    fn subcell_source_expands_and_contracts_without_an_integer_cell_assumption() {
        let grow = super::remap(0.25, 4, 2.0, 1.0, 8.0, 1.0, || false).unwrap();
        close(grow.height_m, 4.0);
        close(grow.source_mass_kg, 2.0);
        close(grow.source_momentum_kg_m_s, 3.5);
        close(grow.source_kinetic_energy_j, 49.0 / 12.0);
        assert_eq!(grow.cells.len(), 2);
        totals(&grow);
        let shrink = super::remap(0.25, 4, 2.0, 1.0, 1.0, 8.0, || false).unwrap();
        assert_eq!(shrink.cells.len(), 1);
        close(shrink.cells[0].fill_fraction, 0.03125);
        close(shrink.source_mass_kg, 0.25);
        totals(&shrink);
        let steady = super::remap(0.25, 4, 2.0, 1.0, 8.0, 8.0, || false).unwrap();
        close(steady.cells[0].fill_fraction, 0.25);
        assert_eq!(steady.source_momentum_kg_m_s, 0.0);
        assert_eq!(steady.source_kinetic_energy_j, 0.0);
    }

    #[test]
    fn successive_partial_density_maps_preserve_closed_mass_and_return_geometry() {
        let mut depth = 2.375;
        let mut cell_mass = 7.0;
        let mass = depth * cell_mass;
        for next_mass in [5.0, 9.0, 3.0, 7.0] {
            let r = super::remap(depth, 16, 0.2, 0.1, cell_mass, next_mass, || false).unwrap();
            close(r.source_mass_kg, mass);
            totals(&r);
            depth = r.cells.iter().map(|c| c.fill_fraction).sum();
            close(depth * next_mass, mass);
            cell_mass = next_mass;
        }
        close(depth, 2.375);
    }

    #[test]
    fn partial_source_keeps_representable_totals_when_density_ratio_overflows() {
        assert!((1e200_f64 / 1e-110).is_infinite());
        let r = super::remap(1e-310, 2, 1e10, 1e10, 1e200, 1e-110, || false).unwrap();
        close(r.source_mass_kg, 1e-110);
        close(r.height_m, 1e10);
        close(r.source_momentum_kg_m_s, 0.5e-110);
        close(r.source_kinetic_energy_j, 1e-110 / 6.0);
        totals(&r);
    }

    #[test]
    fn invalid_partial_depth_is_atomic() {
        for depth in [0.0, -0.25, 2.5, f64::NAN, f64::INFINITY] {
            assert!(super::remap(depth, 2, 1.0, 1.0, 1.0, 1.0, || false).is_none());
        }
        assert!(super::remap(0.5, 0, 1.0, 1.0, 1.0, 1.0, || false).is_none());
        assert!(super::remap(0.5, 2, 1.0, 1.0, 1.0, 1.0, || true).is_none());
    }
}
