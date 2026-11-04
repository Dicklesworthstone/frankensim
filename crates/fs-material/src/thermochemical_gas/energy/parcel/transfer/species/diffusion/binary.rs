//! Exponential integration of the original binary molar-frame diffusion law.
//!
//! For fixed face conductance `G = c D A / L` [mol/s] and parcel mole counts
//! `N_l, N_r`, the composition difference obeys
//! `d(x_l-x_r)/dt = -G (1/N_l + 1/N_r) (x_l-x_r)`.
//! Integrating this equation gives an equal-and-opposite molar exchange that
//! approaches equilibrium without overshooting or a diffusion CFL. Unequal
//! volumes and molecular masses retain their original physical inventories.
//!
//! This is exact for an isothermal, isobaric binary pair with constant source
//! diffusivity. Otherwise it is a frozen-coefficient, first-order step using
//! the same initial face state and upstream partial enthalpies as `diffuse_with`.
//! Source changes and nonisothermal evolution still require time refinement.
//! The original finite-withdrawal, source, mass, energy and entropy gates remain
//! authoritative. No thermal reservoir, extra heat or bulk momentum is added.

use super::{
    BinaryDiffusionMaterial, ClosedGasParcel, GasDiffusionCandidate, GasDiffusionFlux,
    GasDiffusionGeometry, GasDiffusionSolve, GasEnergyError, GasSpeciesExchangeCandidate,
    SpeciesMassTransfer, catalog, finish_candidate, fractions, invalid, poll, positive,
    source_error,
};
use fs_math::det;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Actual integrated binary exchange, distinct from the initial face flux.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryDiffusionStep {
    /// Requested physical interval [s].
    pub duration_s: f64,
    /// `1 - exp(-G (1/N_l + 1/N_r) duration_s)`, in [0,1].
    pub relaxation_fraction: f64,
    /// Integrated transfer of the first canonical species [mol], left to right.
    /// The second species carries exactly its negative before mass conversion.
    pub first_species_transfer_mol: f64,
}

/// Both original thermochemical endpoints and their complete physical account.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryDiffusionCandidate {
    diffusion: GasDiffusionCandidate,
    step: BinaryDiffusionStep,
}

impl BinaryDiffusionCandidate {
    /// Original conservative species and partial-enthalpy exchange.
    pub const fn exchange(&self) -> &GasSpeciesExchangeCandidate {
        self.diffusion.exchange()
    }

    /// Instantaneous flux at the initial face state, including source receipts.
    /// Multiplying this by duration does NOT give the exponential transfer.
    pub const fn initial_flux(&self) -> &GasDiffusionFlux {
        self.diffusion.flux()
    }

    /// Integrated transfer and the actual relaxation fraction.
    pub const fn step(&self) -> &BinaryDiffusionStep {
        &self.step
    }

    /// Actual source-relative entropy change of the complete closed pair [J/K].
    pub const fn entropy_change_j_k(&self) -> f64 {
        self.diffusion.entropy_change_j_k()
    }

    /// Consume the proposal without rebuilding either gas state.
    pub fn into_parts(self) -> (GasDiffusionCandidate, BinaryDiffusionStep) {
        (self.diffusion, self.step)
    }
}

impl ClosedGasParcel {
    /// Advance one binary molar-average diffusion pair by positive relaxation.
    ///
    /// Exactly two species in the joined original catalog and one binary source
    /// are required for an open face. Either species may initially be absent
    /// from one parcel. A closed face retains the original source-free identity.
    /// The canonical source catalog owns species ordering and caloric identity.
    ///
    /// Face T/p, D and the upstream partial enthalpies are sampled at the initial
    /// state. An isothermal/isobaric constant-D pair therefore has the analytic
    /// exponential endpoint, including at intervals that explicit Euler refuses.
    /// Otherwise coefficients are frozen only for this interval; this is not a
    /// fully implicit thermal solve or a time-accuracy certificate. The original
    /// caller T/p limits still apply, and the selected binary source must also
    /// support both final states. All original withdrawal and balance tolerances
    /// are retained, so a caller's small withdrawal allowance can still refuse
    /// a large physically valid exchange.
    ///
    /// Only one shared flux evaluation and one shared two-volume reconstruction
    /// are needed; work does not grow with the diffusion rate times the interval.
    /// Failure or cancellation publishes neither endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_binary_exponential_with(
        &self,
        right: &Self,
        geometry: GasDiffusionGeometry,
        duration_s: f64,
        pairs: &[BinaryDiffusionMaterial<'_>],
        controls: GasDiffusionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryDiffusionCandidate, GasEnergyError> {
        if !positive(duration_s) {
            return Err(invalid(
                "positive finite binary diffusion duration required",
            ));
        }
        let flux = self.diffusion_flux_to(right, geometry, pairs, controls, &mut cancelled)?;
        if geometry.area_m2 == 0.0 {
            return finish_candidate(self, right, flux, &[], controls, &mut cancelled).map(
                |diffusion| BinaryDiffusionCandidate {
                    diffusion,
                    step: BinaryDiffusionStep {
                        duration_s,
                        relaxation_fraction: 0.0,
                        first_species_transfer_mol: 0.0,
                    },
                },
            );
        }
        if flux.species.len() != 2 || flux.pairs.len() != 1 {
            return Err(invalid(
                "exponential diffusion requires exactly two joined gas species",
            ));
        }
        let catalog = catalog(self, right)?;
        let components = catalog.thermodynamics().components();
        let left_moles = self.mass_kg()
            / self
                .state()
                .thermodynamics()
                .properties()
                .molar_mass()
                .value();
        let right_moles = right.mass_kg()
            / right
                .state()
                .thermodynamics()
                .properties()
                .molar_mass()
                .value();
        if !positive(left_moles) || !positive(right_moles) {
            return Err(invalid(
                "binary diffusion needs finite positive molar inventories",
            ));
        }
        // This harmonic capacity is <= each inventory. Never form N_l*N_r or
        // their possibly overflowing sum to obtain the bounded equilibrium.
        let small = left_moles.min(right_moles);
        let large = left_moles.max(right_moles);
        let capacity = small / (1.0 + small / large);
        if !positive(capacity) {
            return Err(invalid("binary diffusion capacity is unrepresentable"));
        }
        let left_x = fractions(self);
        let right_x = fractions(right);
        let first = components[0].species();
        let gap =
            left_x.get(first).copied().unwrap_or(0.0) - right_x.get(first).copied().unwrap_or(0.0);
        let diffusivity = flux.pairs[0].diffusivity_m2_s;
        let mut exponent = (flux.face_pressure_pa / (R * flux.face_temperature_k))
            * diffusivity
            * (geometry.area_m2 / geometry.distance_m)
            * (duration_s / capacity);
        if !exponent.is_normal() {
            // Retain a finite rate-time product when individual scales do not
            // fit; an overflowing positive exponent simply means equilibrium.
            exponent = det::exp(
                det::ln(flux.face_pressure_pa) - det::ln(R) - det::ln(flux.face_temperature_k)
                    + det::ln(diffusivity)
                    + det::ln(geometry.area_m2)
                    - det::ln(geometry.distance_m)
                    + det::ln(duration_s)
                    - det::ln(capacity),
            );
        }
        let relaxation = -det::expm1(-exponent);
        let transfer = (capacity * gap) * relaxation;
        if !relaxation.is_finite()
            || !(0.0..=1.0).contains(&relaxation)
            || !transfer.is_finite()
            || (gap != 0.0 && transfer == 0.0)
        {
            return Err(invalid(
                "binary exponential molar transfer is unrepresentable",
            ));
        }
        let mut rows = Vec::with_capacity(2);
        for (i, species) in components.iter().enumerate() {
            poll(&mut cancelled)?;
            let mass = (if i == 0 { transfer } else { -transfer }) * species.molar_mass().value();
            if !mass.is_finite() || (transfer != 0.0 && mass == 0.0) {
                return Err(invalid(
                    "binary exponential species mass is unrepresentable",
                ));
            }
            rows.push(SpeciesMassTransfer {
                species: species.species().clone(),
                mass_left_to_right_kg: mass,
            });
        }
        let diffusion = finish_candidate(self, right, flux, &rows, controls, &mut cancelled)?;
        // A large stable step must not bypass the source's actual T/p support.
        for parcel in [diffusion.exchange().left(), diffusion.exchange().right()] {
            poll(&mut cancelled)?;
            let state = parcel.state().state();
            let point = crate::thermochemical_gas::energy::updated_point(
                pairs[0].context,
                state.temperature,
                state.pressure,
            )?;
            diffusion.flux().pairs[0]
                .material
                .at_query_point(pairs[0].material, &point)
                .map_err(source_error)?;
        }
        poll(&mut cancelled)?;
        Ok(BinaryDiffusionCandidate {
            diffusion,
            step: BinaryDiffusionStep {
                duration_s,
                relaxation_fraction: relaxation,
                first_species_transfer_mol: transfer,
            },
        })
    }
}
