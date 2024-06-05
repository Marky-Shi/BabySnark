use crate::{common::*, prover::Proof};
use lambdaworks_math::{
    cyclic_group::IsGroup,
    elliptic_curve::traits::{IsEllipticCurve, IsPairing},
    msm::pippenger::msm,
};
use std::ops::Mul;

pub struct VerifyingKey {
    // Ui(τ) * g1, 0 <= i < l
    pub u_tau_g1: Vec<G1Point>,
    // Ui(τ) * g2, 0 <= i < l
    pub u_tau_g2: Vec<G2Point>,
    // t(τ) * g2
    pub t_tau_g2: G2Point,
    // e(g1, g2)^-1
    pub inv_pairing_g1_g2: PairingOutput,
    // β * γ * g1
    pub beta_gamma_g1: G1Point,
    // γ * g2
    pub gamma_g2: G2Point,
}

pub fn verify(proof: &Proof, inputs: &[FrElement], vk: &VerifyingKey) -> bool {
    let v_w = &proof.v_w;
    let v_w_prime = &proof.v_w_prime;
    let h = &proof.h;
    let b_w = &proof.b_w;

    let mut accept = true;

    accept &= Pairing::compute(b_w, &vk.gamma_g2) == Pairing::compute(&vk.beta_gamma_g1, v_w_prime);
    accept &= Pairing::compute(v_w, &TwistCurve::generator())
        == Pairing::compute(&Curve::generator(), v_w_prime);

    let v_u = msm(
        &inputs
            .iter()
            .map(|x| x.representative())
            .collect::<Vec<_>>(),
        &vk.u_tau_g1,
    )
    .unwrap();

    let v_u_prime = msm(
        &inputs
            .iter()
            .map(|elem| elem.representative())
            .collect::<Vec<_>>(),
        &vk.u_tau_g2,
    )
    .unwrap();

    accept &= Pairing::compute(&v_u.operate_with(v_w), &v_u_prime.operate_with(v_w_prime))
        .unwrap()
        .mul(&vk.inv_pairing_g1_g2)
        == Pairing::compute(h, &vk.t_tau_g2).unwrap();

    accept
}
