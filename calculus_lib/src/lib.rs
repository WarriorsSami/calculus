use either::Either;
use num::complex::Complex64;

pub fn newton_secant(x_prev: f32, x_curr: f32, f: fn(f32) -> f32) -> f32 {
    x_curr - f(x_curr) * (x_curr - x_prev) / (f(x_curr) - f(x_prev))
}

pub const EPS: f64 = 1e-3;

pub type Solution = Either<f64, Complex64>;

pub fn solve_linear_equation(a: f64, b: f64) -> Solution {
    Either::Left(-b / a)
}

pub fn solve_quadratic_equation(a: f64, b: f64, c: f64) -> (Solution, Solution) {
    let d = b * b - 4.0 * a * c;

    match d {
        d if d > 0.0 => {
            let x1 = (-b + d.sqrt()) / (2.0 * a);
            let x2 = (-b - d.sqrt()) / (2.0 * a);
            (Either::Left(x1), Either::Left(x2))
        }
        d if d < 0.0 => {
            let x1 = Complex64::new(-b / (2.0 * a), (-d).sqrt() / (2.0 * a));
            let x2 = Complex64::new(-b / (2.0 * a), -(-d).sqrt() / (2.0 * a));
            (Either::Right(x1), Either::Right(x2))
        }
        _ => {
            let x = -b / (2.0 * a);
            (Either::Left(x), Either::Left(x))
        }
    }
}

// p(x) = a[0] * x^n + a[1] * x^(n-1) + ... + a[n-1] * x + a[n]
pub fn solve_polynomial(a: Vec<f64>, x: f64) -> f64 {
    a.iter().enumerate().fold(0.0, |acc, (i, &ai)| {
        acc + ai * x.powi((a.len() - 1 - i) as i32)
    })
}

#[derive(Debug)]
pub struct Matrix {
    data: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn new(data: Vec<Vec<f64>>) -> Self {
        // prepend 0.0 to each row and column
        let mut data = data;
        for row in &mut data {
            row.insert(0, 0.0);
        }
        data.insert(0, vec![0.0; data[0].len()]);

        Self { data }
    }

    pub fn n(&self) -> usize {
        self.data.len() - 1
    }

    pub fn m(&self) -> usize {
        self.data[0].len() - 1
    }

    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.data[i][j]
    }

    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        self.data[i][j] = value;
    }

    pub fn swap_rows(&mut self, i: usize, j: usize) {
        self.data.swap(i, j);
    }

    pub fn swap_columns(&mut self, i: usize, j: usize) {
        for row in &mut self.data {
            row.swap(i, j);
        }
    }

    pub fn get_row_of_max_col_pivot(&self, k: usize) -> usize {
        let mut max = k;
        for i in k + 1..=self.n() {
            if self.at(i, k).abs() > self.at(max, k).abs() {
                max = i;
            }
        }
        max
    }

    pub fn get_row_col_of_max_matrix_pivot(&self, k: usize) -> (usize, usize) {
        let mut max = (k, k);
        for i in k..=self.n() {
            for j in k..=self.m() {
                if self.at(i, j).abs() > self.at(max.0, max.1).abs() {
                    max = (i, j);
                }
            }
        }
        max
    }
}

impl std::fmt::Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for i in 1..=self.n() {
            for j in 1..=self.m() {
                write!(f, "{:.2} ", self.at(i, j))?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}