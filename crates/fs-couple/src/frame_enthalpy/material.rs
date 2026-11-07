//! Original thermal calibrations for the native frame fiber families.
//!
//! Branch arity stays in the type of every source and prepared law. No branch
//! is padded, merged, discarded or reconstructed from total backstress.

use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::plastic::thermal::chaboche::ThermalChabocheSource;
use fs_material::plastic::thermal::tabulated::ThermalTabulatedSource;
use fs_material::plastic::thermal::{FeedbackError, ThermoplasticSource};
use fs_solid::fiber::FiberLaw;

macro_rules! chaboche_families {
    ($($variant:ident: $count:literal),+ $(,)?) => {
        /// Original material/caloric source for a thermal frame cell.
        ///
        /// A frame may mix measured tabulated isotropic/Prager fibers and
        /// cyclic Chaboche fibers. Each cell keeps one original family and
        /// each fiber retains that family's complete accepted history.
        #[derive(Clone, Debug)]
        pub enum FrameThermalSource<'a> {
            Tabulated(ThermalTabulatedSource<'a>),
            $($variant(ThermalChabocheSource<'a, $count>),)+
        }

        impl<'a> FrameThermalSource<'a> {
            pub fn base(&self) -> &ThermoplasticSource<'a> {
                match self {
                    Self::Tabulated(source) => source.base(),
                    $(Self::$variant(source) => source.base(),)+
                }
            }

            pub fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
                self.base().caloric_source()
            }

            /// Prepare once per original cell and feedback proposal. Native
            /// fibers share the owned interval during all equilibrium trials.
            pub(super) fn prepare(
                &self,
                old_temperature_k: f64,
                temperature_k: f64,
                duration_s: f64,
            ) -> Result<FiberLaw, FeedbackError> {
                match self {
                    Self::Tabulated(source) => source
                        .at_uniaxial_temperatures(
                            old_temperature_k,
                            temperature_k,
                            duration_s,
                        )
                        .map(FiberLaw::ThermalTabulated),
                    $(Self::$variant(source) => source
                        .at_uniaxial_temperatures(
                            old_temperature_k,
                            temperature_k,
                            duration_s,
                        )
                        .map(|law| FiberLaw::ThermalChaboche(law.into())),)+
                }
            }
        }

        impl<'a> From<ThermalTabulatedSource<'a>> for FrameThermalSource<'a> {
            fn from(source: ThermalTabulatedSource<'a>) -> Self {
                Self::Tabulated(source)
            }
        }

        $(impl<'a> From<ThermalChabocheSource<'a, $count>> for FrameThermalSource<'a> {
            fn from(source: ThermalChabocheSource<'a, $count>) -> Self {
                Self::$variant(source)
            }
        })+
    };
}

chaboche_families! {
    ChabocheOne: 1,
    ChabocheTwo: 2,
    ChabocheThree: 3,
    ChabocheFour: 4,
    ChabocheFive: 5,
    ChabocheSix: 6,
    ChabocheSeven: 7,
    ChabocheEight: 8,
}
