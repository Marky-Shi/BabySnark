use std::vec;

use crate::{common::*, scs::SquareConstraintSystem};
use lambdaworks_math::polynomial::Polynomial;

#[derive(Debug)]
pub struct SquareSpanProgram {
    pub number_of_public_inputs: usize,
    pub number_of_constraints: usize,
    pub u_polynomials: Vec<Polynomial<FrElement>>,
    pub matrix: Vec<Vec<FrElement>>,
}

impl SquareSpanProgram {
    pub fn calculate_h_coefficients(
        self,
        input: &[FrElement],
        delta: &FrElement,
    ) -> Vec<FrElement> {
        let offset = &ORDER_R_MINUS_1_ROOT_UNITY;
        let p_degree = 2 * self.number_of_constraints;

        let u_evaluated = self.evaluate_scaled_and_acumulated_u(input, p_degree, offset);
        let (t_evaluated, t_inv_evaluated) = self.evaluate_t(p_degree, offset);

        let h_degree = self.number_of_constraints + 1;

        calculate_h_coefficients(
            u_evaluated,
            t_evaluated,
            t_inv_evaluated,
            h_degree,
            offset,
            delta,
        )
    }

    fn evaluate_scaled_and_acumulated_u(
        &self,
        inputs: &[FrElement],
        degree: usize,
        offset: &FrElement,
    ) -> Vec<FrElement> {
        let scaled_and_accumulated = self
            .u_polynomials
            .iter()
            .zip(inputs)
            .map(|(ploy, coeff)| ploy.mul_with_ref(&Polynomial::new_monomial(coeff.clone(), 0)))
            .reduce(|ploy1, ploy2| ploy1 + ploy2)
            .unwrap();

        Polynomial::evaluate_offset_fft(&scaled_and_accumulated, 1, Some(degree), offset).unwrap()
    }

    fn evaluate_t(&self, degree: usize, offset: &FrElement) -> (Vec<FrElement>, Vec<FrElement>) {
        let ploy_new = Polynomial::new_monomial(FrElement::one(), self.number_of_constraints);
        let t_ploy = ploy_new - FrElement::one();

        let t_eval = Polynomial::evaluate_offset_fft(&t_ploy, 1, Some(degree), offset).unwrap();
        let mut t_inv_eval = t_eval.clone();
        FrElement::inplace_batch_inverse(&mut t_inv_eval).unwrap();

        (t_eval, t_inv_eval)
    }

    pub fn number_of_private_inputs(&self) -> usize {
        self.u_polynomials.len() - self.number_of_public_inputs
    }

    pub fn check_vaild(&self, input: &[FrElement]) -> bool {
        for row in &self.matrix {
            let coef = row
                .iter()
                .zip(input)
                .map(|(coef, input)| coef * input)
                .reduce(|a, b| a + b)
                .unwrap();

            if (&coef * &coef).ne(&FrElement::one()) {
                return false;
            }
        }
        true
    }

    pub fn from_scs(scs: SquareConstraintSystem) -> Self {
        let number_of_constraints = scs.clone().number_of_constraints().next_power_of_two();

        let u_ploy = (0..scs.clone().input_size())
            .map(|poly_index| get_u_ploy_from_scs(&scs, poly_index, number_of_constraints))
            .collect();

        let matrix = scs.constraints;
        Self {
            number_of_public_inputs: scs.number_of_public_inputs,
            number_of_constraints,
            u_polynomials: u_ploy,
            matrix,
        }
    }
}

fn calculate_h_coefficients(
    u_evaluated: Vec<FrElement>,
    t_evaluated: Vec<FrElement>,
    t_inv_evaluated: Vec<FrElement>,
    degree: usize,
    offset: &FrElement,
    delta: &FrElement,
) -> Vec<FrElement> {
    let h_eval: Vec<FrElement> = u_evaluated
        .iter()
        .zip(&t_evaluated)
        .zip(&t_inv_evaluated)
        .map(|((ui, ti), ti_inv)| {
            (ui * ui - FrElement::one()) * ti_inv
                + FrElement::from(2) * delta * ui
                + delta * delta * ti
        })
        .collect();

    let mut h_coefficients = Polynomial::interpolate_offset_fft(&h_eval, offset)
        .unwrap()
        .coefficients()
        .to_vec();

    let pad = vec![FrElement::zero(); degree - h_coefficients.len()];

    h_coefficients.extend(pad);
    h_coefficients
}

fn get_u_ploy_from_scs(
    scs: &SquareConstraintSystem,
    poly_index: usize,
    number_of_constraints: usize,
) -> Polynomial<FrElement> {
    let mut u_poly = if poly_index == 0 {
        vec![FrElement::one(); number_of_constraints]
    } else {
        vec![FrElement::zero(); number_of_constraints]
    };

    for (constraint_index, constraint) in scs.constraints.iter().enumerate() {
        u_poly[constraint_index] = constraint[poly_index].clone();
    }

    Polynomial::interpolate_fft::<FrField>(&u_poly).unwrap()
}
