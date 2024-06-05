pub mod common;
pub mod scs;
pub mod ssp;
pub mod utils;

pub mod prover;
pub mod setup;
pub mod verifier;

pub use prover::{Proof, Prover, ProvingKey};
pub use setup::setup;
pub use verifier::{verify, VerifyingKey};
