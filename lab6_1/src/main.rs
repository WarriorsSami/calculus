use calculus_lib::Matrix;

fn main() {
    let mut A = Matrix::new(vec![
        vec![5.0, 2.0, 1.0, 12.0],
        vec![5.0, -6.0, 2.0, -1.0],
        vec![-4.0, 2.0, 1.0, 3.0],
    ]);

    // init identity matrix
    let mut S = Matrix::identity(3);

    for k in 1..A.n() {
        let (p, q) = A.get_row_col_of_matrix_pivot_in_bounds(k, A.n(), A.m() - 1);
        if A.at(p, q) != 0.0 {
            if p != k {
                A.swap_rows(p, k);
            }
            if q != k {
                A.swap_columns(q, k);
                S.swap_columns(q, k);
            }
            for i in k + 1..=A.n() {
                for j in k + 1..=A.m() {
                    A.set(i, j, A.at(i, j) - A.at(i, k) / A.at(k, k) * A.at(k, j));
                }
                A.set(i, k, 0.0);
            }
        } else {
            println!("Cannot solve system Ax = b. A is singular.");
            return;
        }
    }

    println!("{}", A);

    let mut y = vec![0.0; A.n() + 1];
    let mut x = vec![0.0; A.n() + 1];

    y[A.n()] = A.at(A.n(), A.m()) / A.at(A.n(), A.n());
    for i in (1..A.n()).rev() {
        let mut sum = 0.0;
        for j in i + 1..=A.n() {
            sum += A.at(i, j) * y[j];
        }
        y[i] = (A.at(i, A.m()) - sum) / A.at(i, i);
    }

    for i in 1..=A.n() {
        for j in 1..=A.n() {
            x[i] += S.at(i, j) * y[j];
        }
    }

    println!("Solution:");
    for i in 1..=A.n() {
        println!("x[{}] = {}", i, x[i]);
    }
}
