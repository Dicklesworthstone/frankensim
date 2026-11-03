//! Combined physical filament impulse and interval changes. These geometric
//! vortical moments are independent of the near/far partition and core law.

use fs_exec::reduce::Compensated;

use super::FarField;
use super::topology::{circulation_is_closed, visit_physical_segments};
use crate::filament3d::{FilamentWake, Refusal};

/// Checked linear and angular hydrodynamic impulse of the supplied closed
/// physical filament set, divided by uniform fluid density.
///
/// Coordinates and vector components retain the caller's common Cartesian
/// frame. The angular reference is explicit. A nearby integration anchor and
/// compensated partial sums remain private so changing the reference or
/// differencing two snapshots does not subtract large world-origin moments.
/// This is localized-vorticity impulse, not fluid momentum inside an arbitrary
/// finite control volume. It includes the near sheet's closing lifting line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WakeImpulse {
    origin_m: [f64; 3],
    anchor_m: [f64; 3],
    linear: [Compensated; 3],
    angular_at_anchor: [Compensated; 3],
    angular_at_origin: [f64; 3],
}

/// Change of represented wake impulse and its opposite mean wrench equivalent.
///
/// This is a diagnostic contribution, NOT a complete aircraft load model.
/// A force/torque balance must separately account for omitted bound/body
/// vorticity, potential/added-mass impulse, imposed flow, supports, outer
/// boundaries and any pruning/remeshing change. Images are boundary devices
/// and never physical sources here. No force is automatically applied.
///
/// Both endpoint moments use the reported single FIXED spatial origin. In
/// particular, changing the snapshot reference alone produces no exchange.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WakeImpulseExchange {
    /// Common fixed moment reference [m], in the geometry's Cartesian frame.
    pub origin_m: [f64; 3],
    /// Positive elapsed interval [s].
    pub elapsed_s: f64,
    /// Positive uniform density held over the whole interval [kg/m³].
    pub density_kgpm3: f64,
    /// Later minus earlier physical linear wake impulse [N s].
    pub wake_impulse_change_ns: [f64; 3],
    /// Later minus earlier angular wake impulse about `origin_m` [N m s].
    pub wake_angular_impulse_change_nms: [f64; 3],
    /// `-wake_impulse_change_ns / elapsed_s` [N]; only this wake's equivalent.
    pub mean_reaction_force_n: [f64; 3],
    /// `-wake_angular_impulse_change_nms / elapsed_s` [N m]; same fixed origin.
    pub mean_reaction_moment_nm: [f64; 3],
}

impl WakeImpulse {
    /// Fixed origin used by the angular readout [m].
    #[must_use]
    pub fn origin_m(&self) -> [f64; 3] {
        self.origin_m
    }

    /// `I/rho = 1/2 sum Gamma integral x cross dx` [m⁴/s].
    #[must_use]
    pub fn linear_per_density_m4ps(&self) -> [f64; 3] {
        self.linear.map(Compensated::value)
    }

    /// `J/rho = -1/2 sum Gamma integral |x-origin|² dx` [m⁵/s].
    #[must_use]
    pub fn angular_per_density_m5ps(&self) -> [f64; 3] {
        self.angular_at_origin
    }

    /// Re-express the same moment about another fixed origin, without a source
    /// sweep. `J(origin + a) = J(origin) - a cross I`; linear impulse is unchanged.
    ///
    /// # Errors
    /// A non-finite origin or unrepresentable transported angular moment.
    pub fn about_origin(&self, origin_m: [f64; 3]) -> Result<Self, Refusal> {
        admit_origin(origin_m)?;
        let mut result = *self;
        result.origin_m = origin_m;
        let mut angular = result.angular_at_anchor;
        add_cross(
            &mut angular,
            difference(result.anchor_m, origin_m),
            result.linear,
        );
        result.angular_at_origin = checked_values(angular)?;
        Ok(result)
    }

    /// Difference this later snapshot from `before`, using this snapshot's
    /// origin as one fixed reference for BOTH endpoints. Density is uniform
    /// and constant over the positive elapsed interval. This also handles two
    /// snapshots originally queried about different (e.g. moving-CG) origins.
    ///
    /// The result reports the represented wake's change and its negative
    /// average force/moment equivalent. It cannot determine which body,
    /// support, imposed flow or numerical operation supplied that change.
    /// It does not replace surface loads or establish total fluid/body closure.
    ///
    /// # Errors
    /// Invalid density/time, or unrepresentable difference, impulse or rate.
    pub fn exchange_since(
        &self,
        before: &Self,
        elapsed_s: f64,
        density_kgpm3: f64,
    ) -> Result<WakeImpulseExchange, Refusal> {
        if !elapsed_s.is_finite()
            || elapsed_s <= 0.0
            || !density_kgpm3.is_finite()
            || density_kgpm3 <= 0.0
        {
            return Err(refusal(
                "impulse-exchange-invalid",
                "impulse exchange needs positive finite elapsed time and uniform density",
            ));
        }
        let linear_delta = subtract(self.linear, before.linear);
        let mut angular_delta = subtract(self.angular_at_anchor, before.angular_at_anchor);
        // ΔJ_O = ΔJ_local + (anchor_after-anchor_before) × I_before
        //                  + (anchor_after-O) × ΔI.
        // Form the small differences first, before transporting them. This
        // avoids cancellation between two potentially enormous J_O readouts.
        add_cross(
            &mut angular_delta,
            difference(self.anchor_m, before.anchor_m),
            before.linear,
        );
        add_cross(
            &mut angular_delta,
            difference(self.anchor_m, self.origin_m),
            linear_delta,
        );
        let delta_i = checked_values(linear_delta)?;
        let delta_j = checked_values(angular_delta)?;
        let mut change_i = [0.0; 3];
        let mut change_j = [0.0; 3];
        let mut force = [0.0; 3];
        let mut moment = [0.0; 3];
        for i in 0..3 {
            change_i[i] = scaled(delta_i[i], density_kgpm3)?;
            change_j[i] = scaled(delta_j[i], density_kgpm3)?;
            force[i] = rate(change_i[i], elapsed_s)?;
            moment[i] = rate(change_j[i], elapsed_s)?;
        }
        Ok(WakeImpulseExchange {
            origin_m: self.origin_m,
            elapsed_s,
            density_kgpm3,
            wake_impulse_change_ns: change_i,
            wake_angular_impulse_change_nms: change_j,
            mean_reaction_force_n: force,
            mean_reaction_moment_nm: moment,
        })
    }
}

impl FarField {
    /// Measure all three components of linear AND angular impulse of the
    /// supplied physical near/far wake about `origin_m`, divided by density.
    ///
    /// Exact retained far segments and the near sheet's physical shared edges
    /// are the only sources. Cores, age, multipole caches and image systems do
    /// not enter these geometric moments. Nodal circulation must close within
    /// roundoff, so an open pruned/mismatched interface refuses. Removing an
    /// entire closed component cannot be detected from the remaining state;
    /// callers must account for such changes separately from physical loads.
    ///
    /// The conservative near-edge capacity plus retained far-segment count
    /// must not exceed `max_sources`; far-cell count is separately capped by
    /// the same value, including empty cells. Closure uses O(N) scratch and
    /// O(N log N) deterministic indexing; geometric accumulation is linear.
    /// Cancellation is checked before/after work and every 256 source visits.
    /// No input is mutated, and no partial diagnostic is returned on failure.
    ///
    /// # Errors
    /// Invalid origin/geometry, unbalanced circulation, exhausted source cap,
    /// unrepresentable numerical moments, or cancellation.
    pub fn hydrodynamic_impulse_per_density(
        &self,
        near: &FilamentWake,
        origin_m: [f64; 3],
        max_sources: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WakeImpulse, Refusal> {
        checkpoint(&mut cancelled)?;
        admit_origin(origin_m)?;
        admit_source_bound(near, self, max_sources, &mut cancelled)?;
        near.admit_geometry()?;
        if near.downstream_gamma.as_ref().is_some_and(|g| {
            near.rows.is_empty()
                || g.len() != near.cert.n_stations
                || g.iter().any(|v| !v.is_finite())
        }) {
            return Err(refusal(
                "impulse-topology-invalid",
                "invalid near/far interface",
            ));
        }
        let mut lower = [f64::INFINITY; 3];
        let mut upper = [f64::NEG_INFINITY; 3];
        let mut active = false;
        checked_visit(near, self, &mut cancelled, |a, b, gamma| {
            if !gamma.is_finite() || a.iter().chain(&b).any(|v| !v.is_finite()) {
                return Err(refusal(
                    "impulse-geometry-invalid",
                    "physical segments need finite endpoints and circulation",
                ));
            }
            if gamma != 0.0 {
                active = true;
                for i in 0..3 {
                    lower[i] = lower[i].min(a[i]).min(b[i]);
                    upper[i] = upper[i].max(a[i]).max(b[i]);
                }
            }
            Ok(())
        })?;
        if !circulation_is_closed(near, self, || checkpoint(&mut cancelled))? {
            return Err(refusal(
                "impulse-topology-invalid",
                "retained physical sources have unbalanced vortex endpoints",
            ));
        }
        let anchor = if active {
            std::array::from_fn(|i| 0.5 * lower[i] + 0.5 * upper[i])
        } else {
            [0.0; 3]
        };
        let mut linear = [Compensated::zero(); 3];
        let mut angular = [Compensated::zero(); 3];
        checked_visit(near, self, &mut cancelled, |a, b, gamma| {
            if gamma == 0.0 {
                return Ok(());
            }
            let a = difference(a, anchor);
            let b = difference(b, anchor);
            let contribution = cross(a, b).map(|area| gamma * (0.5 * area));
            let centroid = std::array::from_fn(|i| a[i] / 3.0 + b[i] / 3.0);
            let moment = cross(centroid, contribution);
            if a.iter()
                .chain(&b)
                .chain(&contribution)
                .chain(&moment)
                .any(|v| !v.is_finite())
            {
                return Err(unrepresentable());
            }
            // Each edge closes a triangle through the nearby anchor. The
            // radial sides cancel because the physical endpoint flux closes.
            // Triangle angular impulse is its centroid × vector impulse.
            for i in 0..3 {
                linear[i] = linear[i].accumulate(contribution[i]);
                angular[i] = angular[i].accumulate(moment[i]);
            }
            Ok(())
        })?;
        checked_values(linear)?;
        checked_values(angular)?;
        let result = WakeImpulse {
            origin_m,
            anchor_m: anchor,
            linear,
            angular_at_anchor: angular,
            angular_at_origin: [0.0; 3],
        }
        .about_origin(origin_m)?;
        checkpoint(&mut cancelled)?;
        Ok(result)
    }
}

fn admit_source_bound(
    near: &FilamentWake,
    far: &FarField,
    cap: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Refusal> {
    let mut bound = if near.rows.is_empty() {
        0
    } else {
        near.cert
            .n_stations
            .checked_mul(2)
            .and_then(|v| v.checked_add(1))
            .and_then(|edges| near.rows.len().checked_mul(edges))
            .and_then(|v| v.checked_add(near.cert.n_stations))
            .ok_or_else(work_refusal)?
    };
    if far.cells.len() > cap || bound > cap {
        return Err(work_refusal());
    }
    for cell in &far.cells {
        checkpoint(cancelled)?;
        bound = bound
            .checked_add(cell.segs.len())
            .ok_or_else(work_refusal)?;
        if bound > cap {
            return Err(work_refusal());
        }
    }
    Ok(())
}

fn checked_visit(
    near: &FilamentWake,
    far: &FarField,
    cancelled: &mut impl FnMut() -> bool,
    mut visit: impl FnMut([f64; 3], [f64; 3], f64) -> Result<(), Refusal>,
) -> Result<(), Refusal> {
    let mut count = 0_usize;
    let mut result = Ok(());
    visit_physical_segments(near, far, |a, b, gamma| {
        if result.is_err() {
            return;
        }
        if count % 256 == 0 {
            result = checkpoint(cancelled);
        }
        count += 1;
        if result.is_ok() {
            result = visit(a, b, gamma);
        }
    });
    result
}

fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn subtract(a: [Compensated; 3], b: [Compensated; 3]) -> [Compensated; 3] {
    std::array::from_fn(|i| a[i].accumulate(-b[i].sum).accumulate(-b[i].comp))
}

fn add_cross(result: &mut [Compensated; 3], arm: [f64; 3], value: [Compensated; 3]) {
    for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
        result[i] = result[i]
            .accumulate(arm[j] * value[k].sum)
            .accumulate(arm[j] * value[k].comp)
            .accumulate(-arm[k] * value[j].sum)
            .accumulate(-arm[k] * value[j].comp);
    }
}

fn checked_values(values: [Compensated; 3]) -> Result<[f64; 3], Refusal> {
    let result = values.map(Compensated::value);
    if result.iter().any(|v| !v.is_finite()) {
        Err(unrepresentable())
    } else {
        Ok(result)
    }
}

fn scaled(value: f64, density: f64) -> Result<f64, Refusal> {
    let result = value * density;
    if !result.is_finite() || (value != 0.0 && result == 0.0) {
        Err(unrepresentable())
    } else {
        Ok(result)
    }
}

fn rate(impulse: f64, elapsed: f64) -> Result<f64, Refusal> {
    let result = -impulse / elapsed;
    if !result.is_finite() || (impulse != 0.0 && result == 0.0) {
        Err(unrepresentable())
    } else {
        Ok(result)
    }
}

fn admit_origin(origin: [f64; 3]) -> Result<(), Refusal> {
    if origin.iter().any(|v| !v.is_finite()) {
        Err(refusal(
            "impulse-origin-invalid",
            "impulse needs a finite reference origin",
        ))
    } else {
        Ok(())
    }
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(refusal(
            "impulse-cancelled",
            "physical wake impulse query cancelled",
        ))
    } else {
        Ok(())
    }
}

fn work_refusal() -> Refusal {
    refusal(
        "impulse-work-exhausted",
        "physical source capacity exceeds the impulse-query cap",
    )
}

fn unrepresentable() -> Refusal {
    refusal(
        "impulse-unrepresentable",
        "impulse, transport or exchange exceeds the numerical domain",
    )
}

fn refusal(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "use a closed finite physical wake, representable moments and an admitted source cap"
                .into(),
        ],
    }
}
