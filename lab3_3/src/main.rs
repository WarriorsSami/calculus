use calculus_lib::Matrix;

fn main() {
    let mut A = Matrix::new(vec![
        vec![1.0, -1.0, 2.0],
        vec![2.0, 1.0, 1.0],
        vec![1.0, 2.0, -3.0],
    ]);

    for k in 1..A.n() {
        let (p, q) = A.get_row_col_of_max_matrix_pivot(k);
        if A.at(p, q) != 0.0 {
            if p != k {
                A.swap_rows(p, k);
            }
            if q != k {
                A.swap_columns(q, k);
            }
            for i in k + 1..=A.n() {
                for j in k + 1..=A.m() {
                    A.set(i, j, A.at(i, j) - A.at(i, k) / A.at(k, k) * A.at(k, j));
                }
                A.set(i, k, 0.0);
            }
        }
    }

    println!("{}", A);
}

