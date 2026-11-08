//! Original passive structural memory in the resolved fluid/solid equations.
//!
//! A supplied row c [1/s] and time tau [s] retain the energy coordinate
//! x=c.q-z [sqrt(J)]: x_dot=c.v-x/tau, resisting force=c*x,
//! stored energy=x²/2 and heat rate=x²/tau. The original wall modes must use
//! RELAXED stiffness; the memory's unrelaxed stiffness is not added again.
//! For exp(+i omega t), eliminating x gives c*c^T/(1/tau+i omega).
//! The complete matrix is cached once per prepared frequency. Time evolution
//! retains x itself, including its supplied initial history.

use super::super::{PreparedThermoviscous, ThermoviscousError, ThermoviscousFe, finite, poll};
use super::MovingWallField;
use fs_math::{c64::C64, dd::Dd, det};

/// Complete source-channel resources; no channel is removed to satisfy them.
#[derive(Clone, Copy, Debug)]
pub struct MovingWallRelaxationLimits {
    pub maximum_channels: usize,
    /// Retained row/time payload plus one prepared matrix and channel gains.
    /// Counts incoming vector capacities; allocator bookkeeping and native
    /// fluid factorization/time-cache storage and returned field storage are
    /// separate.
    pub maximum_workspace_bytes: usize,
    /// Maximum channel/mode-pair products per attachment or frequency prepare.
    pub maximum_projection_terms: usize,
}

pub(in crate::acoustics::thermoviscous) struct MovingWallRelaxation {
    rows: Vec<Vec<f64>>,
    times: Vec<f64>,
}

#[derive(Default)]
pub(in crate::acoustics::thermoviscous) struct PreparedWallRelaxation {
    pub(in crate::acoustics::thermoviscous) impedance: Vec<C64>,
    gains_s: Vec<C64>,
}

fn limit(what: &'static str, requested: usize, maximum: usize) -> Result<(), ThermoviscousError> {
    if requested > maximum {
        Err(ThermoviscousError::Budget {
            what,
            requested,
            limit: maximum,
        })
    } else {
        Ok(())
    }
}
fn allocation<T>(count: usize) -> Result<Vec<T>, ThermoviscousError> {
    let mut output = Vec::new();
    output.try_reserve_exact(count).map_err(|_| {
        ThermoviscousError::InvalidInput("unable to allocate bounded structural-memory cache")
    })?;
    Ok(output)
}

impl ThermoviscousFe<'_> {
    /// Attach the original complete passive memory in the wall's SAME signed,
    /// mass-normalized basis. Existing instantaneous damping may coexist when
    /// it represents a different loss; these branches must be excluded from it.
    /// Each row supplies c in x_dot=c.v-x/tau. The source owns this law and its
    /// validity band. This method performs no fitting, eigensolve or basis change.
    ///
    /// The structural frequencies must describe relaxed stiffness. Initial x
    /// is supplied explicitly by the causal state; zero_state means relaxed
    /// zero displacement and history, not an inferred preload. A second memory
    /// attachment refuses. An empty attachment explicitly has no memory states.
    pub fn with_moving_wall_relaxation_channels(
        mut self,
        rows_per_s: Vec<Vec<f64>>,
        relaxation_times_s: Vec<f64>,
        limits: MovingWallRelaxationLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousError> {
        poll(cancelled)?;
        let modes = self
            .moving_wall
            .as_ref()
            .ok_or(ThermoviscousError::InvalidInput(
                "structural memory requires the original moving wall",
            ))?
            .structure
            .modes()
            .len();
        let channels = rows_per_s.len();
        if self.moving_wall_relaxation.is_some() || relaxation_times_s.len() != channels {
            return Err(ThermoviscousError::InvalidInput(
                "one complete structural-memory owner and one time per channel required",
            ));
        }
        limit(
            "structural memory channels",
            channels,
            limits.maximum_channels,
        )?;
        let extent = || ThermoviscousError::InvalidInput("structural memory extent overflow");
        let pairs = modes
            .checked_mul(modes)
            .and_then(|n| n.checked_mul(channels))
            .ok_or_else(extent)?;
        limit(
            "structural memory projection",
            pairs,
            limits.maximum_projection_terms,
        )?;
        let mut bytes = rows_per_s
            .capacity()
            .checked_mul(size_of::<Vec<f64>>())
            .and_then(|n| {
                relaxation_times_s
                    .capacity()
                    .checked_mul(size_of::<f64>())
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| {
                channels
                    .checked_mul(size_of::<C64>())
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| {
                modes
                    .checked_mul(modes)
                    .and_then(|m| m.checked_mul(if channels == 0 { 0 } else { size_of::<C64>() }))
                    .and_then(|m| n.checked_add(m))
            })
            .ok_or_else(extent)?;
        for (row, &tau) in rows_per_s.iter().zip(&relaxation_times_s) {
            poll(cancelled)?;
            bytes = row
                .capacity()
                .checked_mul(size_of::<f64>())
                .and_then(|n| bytes.checked_add(n))
                .ok_or_else(extent)?;
            if row.len() != modes
                || row.iter().any(|x| !x.is_finite())
                || !tau.is_finite()
                || tau <= 0.0
                || !(1.0 / tau).is_finite()
            {
                return Err(ThermoviscousError::InvalidInput(
                    "finite complete original memory rows and representable positive relaxation rates required",
                ));
            }
        }
        limit(
            "structural memory workspace",
            bytes,
            limits.maximum_workspace_bytes,
        )?;
        // Admit the original finite unrelaxed stiffness, including every signed
        // cross term. It is checked here but never added to relaxed K separately.
        for i in 0..modes {
            for j in 0..=i {
                let mut value = Dd::ZERO;
                for (k, row) in rows_per_s.iter().enumerate() {
                    if k.is_multiple_of(256) {
                        poll(cancelled)?;
                    }
                    value = value + Dd::from_f64(row[i]) * Dd::from_f64(row[j]);
                }
                let value = value.to_f64();
                if !value.is_finite()
                    || (i == j && value == 0.0 && rows_per_s.iter().any(|row| row[i] != 0.0))
                {
                    return Err(ThermoviscousError::NonFinite);
                }
            }
        }
        poll(cancelled)?;
        self.moving_wall_relaxation = Some(MovingWallRelaxation {
            rows: rows_per_s,
            times: relaxation_times_s,
        });
        Ok(self)
    }

    pub fn moving_wall_relaxation_channels_per_s(&self) -> Option<&[Vec<f64>]> {
        self.moving_wall_relaxation
            .as_ref()
            .map(|x| x.rows.as_slice())
    }
    pub fn moving_wall_relaxation_times_s(&self) -> Option<&[f64]> {
        self.moving_wall_relaxation
            .as_ref()
            .map(|x| x.times.as_slice())
    }

    pub(in crate::acoustics::thermoviscous) fn prepare_wall_relaxation(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedWallRelaxation, ThermoviscousError> {
        poll(cancelled)?;
        let Some(owner) = &self.moving_wall_relaxation else {
            return Ok(PreparedWallRelaxation::default());
        };
        if owner.rows.is_empty() {
            return Ok(PreparedWallRelaxation::default());
        }
        let n = self
            .moving_wall
            .as_ref()
            .expect("admitted moving wall")
            .structure
            .modes()
            .len();
        let mut gains_s = allocation(owner.times.len())?;
        for &tau in &owner.times {
            poll(cancelled)?;
            // Original range-safe complex division; no omega*tau overflow.
            let gain = C64::ONE / C64::new(1.0 / tau, omega);
            if !finite(gain) {
                return Err(ThermoviscousError::NonFinite);
            }
            gains_s.push(gain);
        }
        let mut impedance = allocation(n * n)?;
        impedance.resize(n * n, C64::ZERO);
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..=i {
                let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
                for (k, (row, gain)) in owner.rows.iter().zip(&gains_s).enumerate() {
                    if k.is_multiple_of(256) {
                        poll(cancelled)?;
                    }
                    let product = Dd::from_f64(row[i]) * Dd::from_f64(row[j]);
                    re = re + product * Dd::from_f64(gain.re);
                    im = im + product * Dd::from_f64(gain.im);
                }
                let value = C64::new(re.to_f64(), im.to_f64());
                if !finite(value) {
                    return Err(ThermoviscousError::NonFinite);
                }
                impedance[i * n + j] = value;
                impedance[j * n + i] = value;
            }
        }
        poll(cancelled)?;
        Ok(PreparedWallRelaxation { impedance, gains_s })
    }
}

impl PreparedThermoviscous<'_, '_> {
    pub(in crate::acoustics::thermoviscous) fn wall_relaxation_field(
        &self,
        field: &mut MovingWallField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ThermoviscousError> {
        let Some(owner) = &self.model.moving_wall_relaxation else {
            return Ok(());
        };
        let mut coordinates = allocation(owner.rows.len())?;
        let (mut energy, mut heat) = (Dd::ZERO, Dd::ZERO);
        for ((row, &tau), &gain) in owner
            .rows
            .iter()
            .zip(&owner.times)
            .zip(&self.wall_relaxation.gains_s)
        {
            poll(cancelled)?;
            let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
            for (&c, &v) in row.iter().zip(&field.modal_velocity_m_sqrt_kg_per_s) {
                re = re + Dd::from_f64(c) * Dd::from_f64(v.re);
                im = im + Dd::from_f64(c) * Dd::from_f64(v.im);
            }
            let x = gain * C64::new(re.to_f64(), im.to_f64());
            let dissipative = x.scale(1.0 / det::sqrt(tau));
            let stored = (0.5 * x.re) * (0.5 * x.re) + (0.5 * x.im) * (0.5 * x.im);
            let loss =
                (0.5 * dissipative.re) * dissipative.re + (0.5 * dissipative.im) * dissipative.im;
            if !finite(x) || !stored.is_finite() || !loss.is_finite() {
                return Err(ThermoviscousError::NonFinite);
            }
            coordinates.push(x);
            energy = energy + Dd::from_f64(stored);
            heat = heat + Dd::from_f64(loss);
        }
        field.relaxation_mean_energy_j = energy.to_f64();
        field.relaxation_w = heat.to_f64();
        field.relaxation_coordinates_sqrt_j = coordinates;
        poll(cancelled)?;
        Ok(())
    }
}
